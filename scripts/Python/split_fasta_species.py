######################################################################################
#                                   Split fasta
######################################################################################
# Xinliang Hu
# 2026-04-22
# Python 3.9.6

import re
import os
import logging
from collections import defaultdict
from datetime import datetime

new_dir = "/Users/faa_files/" 
os.chdir(new_dir)

def setup_logger(log_file=None):
    """
    配置日志记录器，同时输出到控制台和文件
    """
    # 创建 logger
    logger = logging.getLogger('fasta_splitter')
    logger.setLevel(logging.INFO)

    # 避免重复添加 handler（多次调用时）
    if logger.handlers:
        return logger

    # 日志格式
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')

    # 控制台输出
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    # 文件输出（如果未指定文件名，自动生成带时间戳的文件名）
    if log_file is None:
        log_file = f"split_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger

def extract_species_core(full_name):
    """
    从方括号内的完整字符串中提取纯物种名（去除菌株编号）
    规则：取前两个空格分隔的词，如果第二个词是 'sp.' 则保留，否则只取前两个词
    """
    parts = full_name.strip().split()
    if len(parts) >= 2:
        genus = parts[0]
        species_epithet = parts[1]
        if species_epithet.lower() == 'sp.':
            return f"{genus} sp"
        else:
            return f"{genus} {species_epithet}"
    return full_name.strip()

def split_fasta_by_species(input_fasta, log_file=None):
    # 初始化日志
    logger = setup_logger(log_file)
    logger.info(f"开始处理文件: {input_fasta}")

    species_entries = defaultdict(list)
    bracket_pattern = re.compile(r'\[([^]]+)\]')

    with open(input_fasta, 'r') as f:
        header = None
        seq_lines = []

        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # 保存前一条记录
                if header is not None:
                    seq = ''.join(seq_lines)
                    match = bracket_pattern.search(header)
                    if match:
                        full_name = match.group(1).strip()
                        species_core = extract_species_core(full_name)
                        species_entries[species_core].append((header, seq))
                    else:
                        logger.warning(f"未找到方括号内容，跳过记录 (行 {line_num}): {header}")
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        # 处理最后一条记录
        if header is not None:
            seq = ''.join(seq_lines)
            match = bracket_pattern.search(header)
            if match:
                full_name = match.group(1).strip()
                species_core = extract_species_core(full_name)
                species_entries[species_core].append((header, seq))
            else:
                logger.warning(f"未找到方括号内容，跳过记录: {header}")

    logger.info(f"共识别出 {len(species_entries)} 个不同的物种")

    # 写入文件
    for species, entries in species_entries.items():
        safe_species = re.sub(r'[\\/*?:"<>| ]', '_', species)
        output_file = f"{safe_species}.fasta"

        with open(output_file, 'w') as out:
            for header, seq in entries:
                out.write(f"{header}\n")
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")

        # 同时输出到屏幕和日志文件
        logger.info(f"已保存: {output_file} (包含 {len(entries)} 条序列)")

    logger.info("处理完成。")

if __name__ == "__main__":
    # 输入文件名
    input_file = "coat_protein.fasta"
    split_fasta_by_species(input_file, log_file="Run.log")
