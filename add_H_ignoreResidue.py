import os
import gzip
import subprocess
import gemmi
import pandas as pd
from tqdm import tqdm
import config

# 设置环境变量
# os.environ['CLIBD_MON'] = 'monomers'
# conda install -c conda-forge gemmi -n xhpi_env


# 目录路径
# root_dir = "/y/peopleconda/bql506/pdb/pdb_mirror/"
# output_dir = "/y/people/bql506/pdb/pdb_mirror_Hadded/"

root_dir = "/Volumes/Sean/pdb_test"
output_dir = "/Volumes/Sean/pdb_test_output"


os.makedirs(output_dir, exist_ok=True)
missing_monomers_file = os.path.join(output_dir, "missing_monomers.csv")

def is_gzip_file(filepath):
    with open(filepath, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

# 用于记录缺失单体信息的DataFrame
missing_monomers_df = pd.DataFrame(columns=["file", "monomer"])

def remove_residue_from_structure(structure, residue_name):
    """
    从结构中移除指定名字的残基。
    
    :param structure: gemmi.Structure 对象
    :param residue_name: 要移除的残基名字
    """
    for model in structure:
        for chain in model:
            residues_to_remove = [i for i, res in enumerate(chain) if res.name == residue_name]
            for i in reversed(residues_to_remove):
                del chain[i]  # 使用 __delitem__ 方法删除残基
    return structure

# 找到根目录中的所有 .gz 文件
gz_files = []
for dirpath, _, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith('.gz') and not filename.startswith('._'):
            gz_files.append(os.path.join(dirpath, filename))

with tqdm(total=len(gz_files), desc="Processing files") as pbar:
    for filepath in gz_files:
        if not is_gzip_file(filepath):
            pbar.update(1)
            continue
        
        try:
            with gzip.open(filepath, 'rb') as file:
                uncompressed_content = file.read().decode('utf-8')
            cif_block = gemmi.cif.read_string(uncompressed_content).sole_block()
            structure = gemmi.make_structure_from_block(cif_block)
        except Exception as e:
            print(f"Error processing CIF content for file {filepath}: {e}")
            pbar.update(1)
            continue

        temp_cif = os.path.join(output_dir, "temp.cif")
        structure.make_mmcif_document().write_file(temp_cif)
        
        temp_output_cif = os.path.join(output_dir, "temp_h.cif")
        command = ['gemmi', 'h', temp_cif, temp_output_cif]
        max_retry = 8  # 最大尝试次数
        retry_count = 0
        success = False
        while retry_count < max_retry:
            retry_count += 1
            result = subprocess.run(command, capture_output=True, text=True)
            
            if result.returncode == 0:
                success = True
                break
            else:
                error_msg = result.stderr.strip()
                if 'Monomer not in the library' in error_msg:
                    monomer_name = error_msg.split('Monomer not in the library: ')[-1].split('.')[0]
                    new_row = pd.DataFrame({"file": [filepath], "monomer": [monomer_name]})
                    missing_monomers_df = pd.concat([missing_monomers_df, new_row], ignore_index=True)
                    try:
                        structure = remove_residue_from_structure(structure, monomer_name)
                        structure.make_mmcif_document().write_file(temp_cif)
                    except Exception as e:
                        print(f"Error processing structure after removing monomer {monomer_name} for file {filepath}: {e}")
                        pbar.update(1)
                        continue
                else:
                    print(f"Failed to add hydrogens for file {filepath}")
                    pbar.update(1)
                    break

        if success:
            try:
                modified_structure = gemmi.read_structure(temp_output_cif)
                cif_string = modified_structure.make_mmcif_document().as_string()
                relative_path = os.path.relpath(filepath, root_dir)
                output_cif_gz = os.path.join(output_dir, relative_path)
                output_cif_gz_dir = os.path.dirname(output_cif_gz)
                if not os.path.exists(output_cif_gz_dir):
                    os.makedirs(output_cif_gz_dir)
                with open(temp_output_cif, 'rt') as temp_output_file:
                    with gzip.open(output_cif_gz, 'wt') as gz_output:
                        gz_output.writelines(temp_output_file)
            except Exception as e:
                print(f"Error writing gzipped CIF to file {output_cif_gz}: {e}")
                pbar.update(1)
        else:
            print(f"Failed to add hydrogens for file {filepath} after {max_retry} retries")
            pbar.update(1)

        os.remove(temp_cif)
        os.remove(temp_output_cif)
        pbar.update(1)
        missing_monomers_df.to_csv(missing_monomers_file, index=False)
