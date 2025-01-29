import pandas as pd
import os
import numpy as np


class Datafile_to_grids():
    def __init__(self,
                 GP_data,      # data file 原始数据
                 grid_size,    # 网格划分的尺寸
                 limits,       # 实际的 x/y/z 边界上下限 
                 limits_new,  # 需要补全到的 x/y/z 边界上下限
                 output_filename = "CHECK_DataFile"
                 ):
        self.output_filename = output_filename
        self.GP_data = GP_data
        self.grid_size = grid_size
        self.xmin, self.xmax = limits["x"]
        self.ymin, self.ymax = limits["y"]
        self.zmin, self.zmax = limits["z"]
        self.xmin_new, self.xmax_new = limits_new["x"]
        self.ymin_new, self.ymax_new = limits_new["y"]
        self.zmin_new, self.zmax_new = limits_new["z"]

        # 计算原始边界和新边界的长度
        self.xlength = self.xmax - self.xmin
        self.ylength = self.ymax - self.ymin
        self.zlength = self.zmax - self.zmin
        self.xlength_new = self.xmax_new - self.xmin_new
        self.ylength_new = self.ymax_new - self.ymin_new
        self.zlength_new = self.zmax_new - self.zmin_new

        # 自动执行
        self.GP_data_pretreat()
        self.GP_data_analysis()
        self.GP_data_aligning()
        self.GP_data_to_grid()
        self.check_atom_data_list()

    def GP_data_pretreat(self):
        self.GP_data = np.delete(self.GP_data, 1, axis=1)  # 删除 atom type 列

    def GP_data_analysis(self):
        self.GP_atom_count = self.GP_data.shape[0]

    def GP_data_aligning(self):
        '''
        对data file 进行补全
        '''
        # 计算需要复制的倍数，包括正方向和负方向
        x_repeats_pos = int(np.ceil((self.xmax_new - self.xmax) / self.xlength))
        x_repeats_neg = int(np.ceil((self.xmin - self.xmin_new) / self.xlength))
        y_repeats_pos = int(np.ceil((self.ymax_new - self.ymax) / self.ylength))
        y_repeats_neg = int(np.ceil((self.ymin - self.ymin_new) / self.ylength))
        z_repeats_pos = int(np.ceil((self.zmax_new - self.zmax) / self.zlength))
        z_repeats_neg = int(np.ceil((self.zmin - self.zmin_new) / self.zlength))

        # 扩展原子坐标，包括正方向和负方向
        extended_atoms = []
        new_atom_id = 1
        for i in range(-x_repeats_neg, x_repeats_pos + 1):
            for j in range(-y_repeats_neg, y_repeats_pos + 1):
                for k in range(-z_repeats_neg, z_repeats_pos + 1):
                    for atom in self.GP_data:
                        new_atom = atom.copy()
                        new_atom[0] = new_atom_id  # 更新原子ID
                        new_atom[1] = atom[1] + i * self.xlength
                        new_atom[2] = atom[2] + j * self.ylength
                        new_atom[3] = atom[3] + k * self.zlength
                        extended_atoms.append(new_atom)
                        new_atom_id += 1

        # 转换为numpy数组
        extended_atoms = np.array(extended_atoms)

        # 过滤掉超出新边界的原子
        extended_atoms = extended_atoms[
            (extended_atoms[:, 1] >= self.xmin_new) & (extended_atoms[:, 1] < self.xmax_new) &
            (extended_atoms[:, 2] >= self.ymin_new) & (extended_atoms[:, 2] < self.ymax_new) &
            (extended_atoms[:, 3] >= self.zmin_new) & (extended_atoms[:, 3] < self.zmax_new)
        ]

        # 更新GP_data，并重新分配ID以确保连续
        self.GP_data = extended_atoms
        self.GP_data[:, 0] = np.arange(1, self.GP_data.shape[0] + 1)
        self.GP_atom_count_extended = self.GP_data.shape[0]

    def GP_data_to_grid(self):
        '''
        将原子数据转换为网格数据
        '''
        # N 格子, N+1 个点
        x_bins = np.linspace(self.xmin_new, self.xmax_new, self.grid_size[0] + 1)
        y_bins = np.linspace(self.ymin_new, self.ymax_new, self.grid_size[1] + 1)
        z_bins = np.linspace(self.zmin_new, self.zmax_new, self.grid_size[2] + 1)

        # 生成 X * Y * Z 形状的矩阵 (grid_size=(XY,Z))
        grid_counts = np.zeros(self.grid_size)

        for atom in self.GP_data:
            x_idx = np.digitize(atom[1], x_bins) - 1    # digitize() 从1开始
            y_idx = np.digitize(atom[2], y_bins) - 1
            z_idx = np.digitize(atom[3], z_bins) - 1
            if x_idx < self.grid_size[0] and y_idx < self.grid_size[1] and z_idx < self.grid_size[2]:
                grid_counts[x_idx, y_idx, z_idx] += 1

        self.grid_counts = grid_counts

        # 打印结果以检查
        #print(f"Grid counts shape: {self.grid_counts.shape}")
        #print(self.grid_counts)

    def check_atom_data_list(self):
        """
        将atom_data_list的内容输出为LAMMPS dump文件格式
        :param output_filename: 输出文件名
        """
        
        with open(self.output_filename, 'w') as data_file:
            # 写入基本信息
            data_file.write("LAMMPS data file via write_data, version 2 Aug 2023, timestep = 50000, units = metal\n\n")
            data_file.write(f"{self.GP_atom_count_extended} atoms\n")
            data_file.write("1 atom types\n\n")
            data_file.write(f"{self.xmin_new} {self.xmax_new} xlo xhi\n") 
            data_file.write(f"{self.ymin_new} {self.ymax_new} ylo yhi\n")
            data_file.write(f"{self.zmin_new} {self.zmax_new} zlo zhi\n\n") 
            data_file.write("Masses\n\n")
            data_file.write("1 12.0107\n\n")
            data_file.write("Atoms # full\n\n")
            for atom in self.GP_data:
                atom_inserted = np.insert(atom, 1, 1)                            # 在第二项位置插入 atom type   为 1
                atom_inserted = np.insert(atom_inserted, 1, 1)                   # 在第二项位置插入 molecule ID 为 1
                atom_inserted = np.insert(atom_inserted, 3, 0)                   # 在第二项位置插入 charge      为 0
                atom_inserted = np.insert(atom_inserted, len(atom_inserted), 0)  # 在第二项位置插入 velocity X  为 0
                atom_inserted = np.insert(atom_inserted, len(atom_inserted), 0)  # 在第二项位置插入 velocity Y  为 0
                atom_inserted = np.insert(atom_inserted, len(atom_inserted), 0)  # 在第二项位置插入 velocity Z  为 0
           
                # 将 atom ID 和 molecule ID 转换为整数, 其他保持原有精度
                formatted_line = f"{int(atom_inserted[0])} {int(atom_inserted[1])} {int(atom_inserted[2])} {int(atom_inserted[3])} " \
                                 f"{atom_inserted[4]} {atom_inserted[5]} {atom_inserted[6]} " \
                                 f"{int(atom_inserted[7])} {int(atom_inserted[8])} {int(atom_inserted[9])}"

                data_file.write(formatted_line + "\n")