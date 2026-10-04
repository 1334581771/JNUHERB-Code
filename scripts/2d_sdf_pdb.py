#!/usr/bin/env python3
import os
import subprocess
import re
import argparse
import time
import sys
from multiprocessing import Pool, cpu_count

# --- 核心依赖检查 ---
try:
    from openbabel import pybel
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("❌ 缺失依赖！请执行: conda install -c conda-forge openbabel rdkit -y")
    sys.exit(1)

def sanitize_filename(name):
    """清理文件名中的非法字符"""
    return re.sub(r'[\\/*?:"<>|]', "_", str(name)).strip()

def run_rdkit_fallback(sdf_block, output_path):
    """RDKit 备选方案：使用 ETKDGv3 算法并进行 UFF 优化"""
    try:
        mol = Chem.MolFromMolBlock(sdf_block, removeHs=False)
        if mol is None: return False
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1: return False
        AllChem.UFFOptimizeMolecule(mol) 
        Chem.MolToPDBFile(mol, output_path)
        return True
    except:
        return False

def worker_task(data):
    """双引擎并行任务：优先 OBabel (--best 精度)，失败则 RDKit"""
    sdf_block, file_id, output_path, ph = data
    
    # 策略 1: OpenBabel
    try:
        cmd = [
            "obabel", "-isdf", "-opdb", 
            "-O", output_path, 
            "--gen3d", "best",      # 借鉴你的代码：高质量构象
            "-p", str(ph),          # 生理 pH 质子化
            "--minimize", "--ff", "MMFF94" # 强制力场优化
        ]
        # 增加 timeout 防止某些复杂大分子导致进程死锁
        subprocess.run(cmd, input=sdf_block.encode(), capture_output=True, check=True, timeout=120)
        if os.path.exists(output_path):
            return True, file_id, "OB"
    except Exception:
        pass

    # 策略 2: RDKit 备选
    if run_rdkit_fallback(sdf_block, output_path):
        return True, file_id, "RDKit"
    
    return False, file_id, "FAIL"

def main():
    # --- 完善参数解析 (借鉴别人的简洁性，保留你的灵活性) ---
    parser = argparse.ArgumentParser(description="分子库预处理流水线：2D SDF 转 3D PDB")
    parser.add_argument("-i", "--input", required=True, help="输入的 2D SDF 文件路径")
    parser.add_argument("-o", "--output", default="prepared_ligands", help="输出 PDB 目录 (默认: prepared_ligands)")
    parser.add_argument("-ph", "--ph", type=float, default=7.4, help="质子化 pH 值 (默认: 7.4)")
    parser.add_argument("-w", "--workers", type=int, default=int(cpu_count() * 0.8), 
                        help="并行核心数 (默认: 总核数的 80%%)")
    args = parser.parse_args()

    if not os.path.exists(args.output): 
        os.makedirs(args.output)

    tasks = []
    print(f"📦 正在解析 {args.input} ...")
    
    try:
        for i, mol in enumerate(pybel.readfile("sdf", args.input)):
            # --- 强化 ID 提取逻辑 (结合两份代码的优点) ---
            fid = None
            # 遍历可能的 ID 字段
            for key in ["CAS","Name","Catalog_NO", "ID", "IDNUMBER"]:
                if key in mol.data:
                    fid = mol.data[key]
                    break
            
            if not fid:
                fid = mol.title.strip() or f"mol_{i+1:04d}"
            
            clean_fid = sanitize_filename(fid)
            out = os.path.join(args.output, f"{clean_fid}.pdb")
            
            # --- 断点续传逻辑 (借鉴别人的代码，避免重复跑相同任务) ---
            if not os.path.exists(out):
                tasks.append((mol.write("sdf"), fid, out, args.ph))
                
    except Exception as e:
        print(f"❌ 读取文件失败: {e}")
        sys.exit(1)

    if not tasks:
        print("✅ 没有需要处理的新任务（所有分子已存在于输出目录中）。")
        return

    print(f"🚀 启动转换: {len(tasks)} 个分子 | 核心数: {args.workers} | pH: {args.ph}")
    
    start_time = time.time()
    suc = 0
    rk_rescue = 0
    fail_list = []

    with Pool(processes=args.workers) as pool:
        # imap_unordered 效率最高
        for idx, (is_ok, fid, engine) in enumerate(pool.imap_unordered(worker_task, tasks)):
            if is_ok:
                suc += 1
                if engine == "RDKit": rk_rescue += 1
            else:
                fail_list.append(fid)
            
            # 定时输出进度
            if (idx + 1) % 50 == 0 or (idx + 1) == len(tasks):
                elapsed = time.time() - start_time
                avg_speed = elapsed / (idx + 1)
                rem_time = avg_speed * (len(tasks) - (idx + 1))
                print(f"进度: {idx+1}/{len(tasks)} | 成功: {suc} (RDKit 补救: {rk_rescue}) | 剩余预估: {rem_time:.1f}s")

    print(f"\n✨ 任务完成！")
    print(f"统计: 成功 {suc} (OpenBabel {suc - rk_rescue}, RDKit {rk_rescue}) | 失败 {len(fail_list)}")
    
    if fail_list:
        log_file = os.path.join(args.output, "conversion_failed.log")
        with open(log_file, "w") as f:
            f.write("\n".join(fail_list))
        print(f"⚠️ 失败列表已记录至: {log_file}")

if __name__ == "__main__":
    main()