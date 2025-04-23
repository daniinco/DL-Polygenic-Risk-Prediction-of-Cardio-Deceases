import os
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path


def extract_pgs_snps_to_dataset(pgs_file_path, plink_file_path, output_dir):
    """
    Извлекает SNP, используемые в PGS и создает датасет
    
    pgs_file_path (str): Путь до файла с PGS
    plink_file_path (str): Путь до PLINK
    output_dir (str): Путь до папки для результатов
    
    Возвращает:
    tuple: (X, y) где X - DataFrame с генотипами SNP, y - Series с фенотипами
    """
    os.makedirs(output_dir, exist_ok=True)
    
    plink_basename = os.path.basename(plink_file_path)
    pgs_basename = os.path.basename(pgs_file_path).replace('.txt', '')
    extracted_plink = os.path.join(output_dir, f"{plink_basename}_{pgs_basename}_extracted")
    
    csv_output = os.path.join(output_dir, f"{plink_basename}_{pgs_basename}_dataset.csv")
    print(f"Обработка PGS: {pgs_file_path}")
    print(f"исходный PLINK: {plink_file_path}")
    
    try:
        pgs_data = pd.read_csv(pgs_file_path, sep='\t', comment='#', 
                             dtype={'rsID': str, 'chr_name': str, 'chr_position': str})
        
        print(f"Загруже PGS с {len(pgs_data)} SNP")
        
        is_rsid_based = 'rsID' in pgs_data.columns
        is_position_based = 'chr_name' in pgs_data.columns and 'chr_position' in pgs_data.columns
        
        if not is_rsid_based and not is_position_based:
            print("В файле нет rsID и позиционной инфорамции")
            return None, None
        
        bim_file = f"{plink_file_path}.bim"
        bim_data = pd.read_csv(bim_file, sep='\s+', header=None, 
                              names=['chr', 'rsID', 'cm', 'pos', 'A1', 'A2'])
        
        pos_dict = {}
        for _, row in bim_data.iterrows():
            pos_id = f"chr{row['chr']}:{row['pos']}"
            pos_dict[pos_id] = row['rsID']
        
        if not is_rsid_based and is_position_based:
            pgs_data['rsID'] = pgs_data.apply(
                lambda row: f"chr{row['chr_name']}:{row['chr_position']}", axis=1)
            
            pgs_data['original_id'] = pgs_data['rsID'].copy()
            pgs_data['rsID'] = pgs_data['rsID'].map(lambda x: pos_dict.get(x, x))
        
        pgs_snps = set(pgs_data['rsID'])
        bim_snps = set(bim_data['rsID'])
        common_snps = pgs_snps.intersection(bim_snps)
        
        print(f"Количество SNP в PGS: {len(pgs_snps)}")
        print(f"Количество общих SNP: {len(common_snps)}")
        print(f"Процент покрытия: {len(common_snps)/len(pgs_snps)*100:.2f}%")
        
        snp_list_file = os.path.join(output_dir, "snp_list.txt")
        with open(snp_list_file, 'w') as f:
            for snp in common_snps:
                f.write(f"{snp}\n")
        
        extract_command = [
            'plink',
            '--bfile', plink_file_path,
            '--extract', snp_list_file,
            '--make-bed',
            '--out', extracted_plink
        ]
        
        extract_result = subprocess.run(extract_command, capture_output=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if extract_result.returncode != 0:
            print(f"Ошикба при извлечении SNP {extract_result.stderr}")
            return None, None
        
        print(f"Создан новый PLINK файл с {len(common_snps)} SNP: {extracted_plink}")
        
        # Используем PLINK --recode A для получения аддитивного кодирования генотипов
        recode_command = [
            'plink',
            '--bfile', extracted_plink,
            '--recode', 'A',
            '--out', extracted_plink
        ]
        
        recode_result = subprocess.run(recode_command, capture_output=True, text=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        
        if recode_result.returncode != 0:
            print(f"Ошибка при перекодировании файла {recode_result.stderr}")
            return None, None
        
        raw_file = f"{extracted_plink}.raw"
        genotypes = pd.read_csv(raw_file, sep='\s+')
        
        print(f"Загружена матрица генотипов размером {genotypes.shape}")
        
        fam_file = f"{plink_file_path}.fam"
        fam_data = pd.read_csv(fam_file, sep='\s+', header=None, 
                              names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE'])
        
        X = genotypes.drop(['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis=1)
        
        genotypes['ID'] = genotypes['FID'].astype(str) + '_' + genotypes['IID'].astype(str)
        fam_data['ID'] = fam_data['FID'].astype(str) + '_' + fam_data['IID'].astype(str)
        
        phenotype_dict = dict(zip(fam_data['ID'], fam_data['PHENOTYPE']))
        y = genotypes['ID'].map(phenotype_dict)

        snp_info = bim_data[bim_data['rsID'].isin(common_snps)].copy()
        snp_info.to_csv(os.path.join(output_dir, f"{plink_basename}_{pgs_basename}_snp_info.csv"), index=False)

        X.to_csv(os.path.join(output_dir, f"{plink_basename}_{pgs_basename}_X.csv"), index=False)
        y.to_csv(os.path.join(output_dir, f"{plink_basename}_{pgs_basename}_y.csv"), index=False)

        full_dataset = pd.concat([genotypes[['FID', 'IID']], X, y.rename('phenotype')], axis=1)
        full_dataset.to_csv(csv_output, index=False)
        
        print(f"Датасет в {csv_output}")
        print(f"X содержит {X.shape[1]} SNP и {X.shape[0]} образцов")
        
        return X, y
        
    except Exception as e:
        print(f"Ошибка {e}")
        return None, None


if __name__ == "__main__":
    print("для извлечения SNP из PLINK файлов на основе PGS")
    print("используйте функцию extract_pgs_snps_to_dataset")
