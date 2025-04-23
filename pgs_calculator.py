import os
import subprocess
import pandas as pd
import numpy as np
import glob
from pathlib import Path


def calculate_pgs(plink_file_path, pgs_folder_path, output_csv_path):
    """
    Рассчитывает полигенные риск скоры для заданного PLINK файла на основе файлов весв
    
    plink_file_path (str): путь до PLINK файла (без расширения)
    pgs_folder_path (str): путь до папки с .txt файлами для расчета PGS
    output_csv_path (str): путь для сохранения результирующего .csv
    
    Возвращает:
    tuple: (X, y) X - DataFrame с PGS, y - Series с фенотипами
    """
    
    temp_dir = os.path.join(os.path.dirname(output_csv_path), 'temp_pgs')
    print(1)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
        print(f"Создана директория для временных файлов: {temp_dir}")
    plink_base_name = os.path.basename(plink_file_path)
    
    bim_file = f"{plink_file_path}.bim"
    bim_data = pd.read_csv(bim_file, sep='\s+', header=None, 
                           names=['chr', 'rsID', 'cm', 'pos', 'A1', 'A2'])
    
    # Создаем множество SNP и словарь для поиска по позиции
    snps = set(bim_data['rsID'])
    pos_dict = {}
    for _, row in bim_data.iterrows():
        pos_id = f"chr{row['chr']}:{row['pos']}"
        pos_dict[pos_id] = row['rsID']
    
    print(f"файл {bim_file}, содержит {len(snps)} SNP")
    pgs_files = glob.glob(os.path.join(pgs_folder_path, '*.txt'))
    print(f"Найдено {len(pgs_files)} PGS")
    
    profile_files = []
    
    # Обработка каждого PGS файла
    for pgs_file_path in pgs_files:
        pgs_name = os.path.basename(pgs_file_path).replace('.txt', '')
        print(f"\n{'='*42}\nОбработка PGS: {pgs_name}")
        
        try:
            pgs_data = pd.read_csv(pgs_file_path, 
                                  sep='\t', 
                                  comment='#', 
                                  dtype={'rsID': str, 'chr_name': str, 'chr_position': str})
            
            print(f"PGS {pgs_name}, содержит {len(pgs_data)} SNP")
            
            # Определяем тип файла (rsID-based или position-based)
            is_rsid_based = 'rsID' in pgs_data.columns
            is_position_based = 'chr_name' in pgs_data.columns and 'chr_position' in pgs_data.columns
            
            required_base_columns = ['effect_allele', 'effect_weight']
            if not is_rsid_based and not is_position_based:
                print("В файле нет rsID и позиционной инфорамции")
                continue
                
            missing_columns = [col for col in required_base_columns if col not in pgs_data.columns]
            if missing_columns:
                print(f"отсутствуют необходимы колонки: {', '.join(missing_columns)}")
                continue
            
            if not is_rsid_based and is_position_based:
                print("Файл содержит позиционную информацию вместо rsID. Создаем синтетиические идентификаторы.")
                pgs_data['rsID'] = pgs_data.apply(lambda row: f"chr{row['chr_name']}:{row['chr_position']}", axis=1)
                
                # Заменяем позиционные идентификаторы на rsID если возможно
                pgs_data['original_id'] = pgs_data['rsID'].copy()
                pgs_data['rsID'] = pgs_data['rsID'].map(lambda x: pos_dict.get(x, x))
            
            common_snps = set(pgs_data['rsID']).intersection(snps)
            
            print(f"Количество SNP в PGS: {len(pgs_data)}")
            print(f"Количество общих SNP: {len(common_snps)}")
            print(f"Процент покрытия: {len(common_snps)/len(pgs_data)*100:.2f}%")
            
            if len(common_snps) == 0:
                print(f"нет общих SNP для {pgs_name}")
                continue
            
            # Создание файла весов для PLINK
            filtered_pgs = pgs_data[pgs_data['rsID'].isin(common_snps)]
            weight_file = os.path.join(temp_dir, f'{pgs_name}_weights.txt')
            
            with open(weight_file, 'w') as f:
                f.write("SNP A1 WEIGHT\n")
                for _, row in filtered_pgs.iterrows():
                    snp_id = row['rsID']
                    effect_allele = row['effect_allele']
                    weight = row['effect_weight']
                    f.write(f"{snp_id} {effect_allele} {weight}\n")
            
            # Расчет PGS
            output_prefix = os.path.join(temp_dir, pgs_name)
            pgs_command = [
                'plink',
                '--bfile', plink_file_path,
                '--score', weight_file, '1', '2', '3', 'header',
                '--out', output_prefix
            ]
            
            print(f"Выполняем расчет PGS {pgs_name}")
            result = subprocess.run(pgs_command, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"PGS {pgs_name} успешно рассчитан")
                profile_files.append(f'{output_prefix}.profile')
            else:
                print(f"Ошибк при расчете PGS {pgs_name}", result.stderr)
                continue
                
        except Exception as e:
            print(f"ошибка при обработке {pgs_name} {e}")
            continue
    
    if not profile_files:
        print("НЕ удалось рассчитать ни один PGS")
        return None, None
    
    # Объединение всех PGS в один датасет
    try:
        base_df = pd.read_csv(profile_files[0], sep='\s+')
        
        combined = base_df[['FID', 'IID']]
        
        for profile_file in profile_files:
            pgs_name = os.path.basename(profile_file)
            pgs_name = pgs_name.replace('.profile', '')
            
            df = pd.read_csv(profile_file, sep='\s+')
            combined[pgs_name] = df['SCORE']
        
        print(f"Объеденены PGS для {len(combined)} образцов и {len(combined.columns)-2} PGS")
        
        fam_file = f"{plink_file_path}.fam"
        fam_data = pd.read_csv(fam_file, sep='\s+', header=None, 
                              names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE'])
        
        phenotype_data = fam_data[['FID', 'IID', 'PHENOTYPE']]
        combined = pd.merge(combined, phenotype_data, on=['FID', 'IID'], how='left')
        
        X = combined.drop(['FID', 'IID', 'PHENOTYPE'], axis=1)
        y = combined['PHENOTYPE']
        
        final_dataset = pd.concat([combined[['FID', 'IID']], X, pd.Series(y, name='y')], axis=1)
        
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)
        final_dataset.to_csv(output_csv_path, index=False)
        print(f"Сохраняю в {output_csv_path}")
        
        for file in glob.glob(os.path.join(temp_dir, '*')):
            try:
                os.remove(file)
            except:
                pass
        
        return X, y
        
    except Exception as e:
        print(f"Ошибка при объединении {e}")
        return None, None


if __name__ == "__main__":
    print("модуль для расчета полигенных рискскоров")
    print("используйте функцию calculate_pgs для расчета")
