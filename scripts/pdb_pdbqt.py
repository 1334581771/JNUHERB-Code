#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse
import signal
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- 超时控制装饰器 ---
def timeout_handler(signum, frame):
    raise TimeoutError("分子处理超时")

# --- 依赖加载 ---
try:
    from openbabel import pybel
    from rdkit import Chem
    from rdkit import RDLogger
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    RDLogger.DisableLog('rdApp.*') 
except ImportError:
    pass

def process_one_file(pdb_path, output_dir, ph_value, timeout_limit=60):
    """
    带超时保护和 API 优化的转换函数
    """
    pdb_path = Path(pdb_path)
    out_path = Path(output_dir) / f"{pdb_path.stem}.pdbqt"
    
    # 1. 断点续传
    if out_path.exists() and out_path.stat().st_size > 0:
        return "SKIP", pdb_path.name

    # 开启超时信号 (仅限 Linux)
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(timeout_limit)

    try:
        # 2. 策略 A: Meeko + RDKit
        try:
            mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
            if mol:
                mol.UpdatePropertyCache(strict=False)
                mol = Chem.AddHs(mol, addCoords=True)
                
                preparer = MoleculePreparation()
                mol_setups = preparer.prepare(mol)
                
                # 检查旋转键数量，防止 Vina 崩溃或 Meeko 无限计算
                if len(preparer.get_active_tours()) > 32:
                    signal.alarm(0) # 清除时钟
                    return "FAIL", f"{pdb_path.name} (Torsions > 32)"

                # 使用新版 Meeko API 写入文件
                if mol_setups:
                    pdbqt_str = PDBQTWriterLegacy.write_string(mol_setups[0])
                    with open(out_path, "w") as f:
                        f.write(pdbqt_str)
                    signal.alarm(0)
                    return "RDKIT", pdb_path.name
        except Exception:
            pass

        # 3. 策略 B: OpenBabel 救火
        try:
            mol = next(pybel.readfile("pdb", str(pdb_path)))
            mol.OBMol.AddHydrogens(False, True, ph_value)
            # 注意：pybel 的 write 会自动处理文件写入
            mol.write("pdbqt", str(out_path), overwrite=True)
            signal.alarm(0)
            return "OB", pdb_path.name
        except Exception:
            pass

    except TimeoutError:
        return "FAIL", f"{pdb_path.name} (Timeout)"
    except Exception as e:
        return "FAIL", f"{pdb_path.name} ({str(e)})"
    finally:
        signal.alarm(0) # 确保闹钟关闭

    return "FAIL", pdb_path.name

def main():
    parser = argparse.ArgumentParser(description="优化版 PDB to PDBQT 处理器")
    parser.add_argument("-i", "--input", default="./pdb", help="PDB输入目录")
    parser.add_argument("-o", "--output", default="./ligand", help="PDBQT输出目录")
    parser.add_argument("-p", "--ph", type=float, default=7.4, help="pH值")
    # 建议 192 核服务器使用 64-128，给 I/O 留余量
    parser.add_argument("-w", "--workers", type=int, default=64)
    args = parser.parse_args()

    in_dir = Path(args.input)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_files = list(in_dir.glob("*.pdb"))
    if not pdb_files:
        print(f"❌ 未找到 PDB 文件")
        return

    print(f"🚀 启动 | 核心: {args.workers} | 待处理: {len(pdb_files)}")

    stats = {"OB": 0, "RDKIT": 0, "FAIL": 0, "SKIP": 0}
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_one_file, f, args.output, args.ph): f for f in pdb_files}
        
        for i, future in enumerate(as_completed(futures), 1):
            try:
                res_type, filename = future.result()
                stats[res_type] += 1
            except Exception as e:
                stats["FAIL"] += 1
            
            if i % 100 == 0 or i == len(pdb_files):
                print(f"进度: {i}/{len(pdb_files)} | 成功: {stats['OB']+stats['RDKIT']} | 跳过: {stats['SKIP']} | 失败: {stats['FAIL']}")

    print(f"\n📊 统计: RDKit:{stats['RDKIT']} | OB:{stats['OB']} | 跳过:{stats['SKIP']} | 失败:{stats['FAIL']}")

if __name__ == "__main__":
    main()