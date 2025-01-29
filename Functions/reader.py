import pandas as pd
import os
import numpy as np

# fetch data from EXCEL
def fetch(len_i=10, 
          sigma_i=14, 
          temp_i=373, 
          item_i='flux\n(L/m^2/h)',
          original_file_path='../../Database/LmpGP.xlsx'
          ):
    # 加载Excel文件
    df = pd.read_excel(original_file_path)
    # 查找行号
    row_number = df[
        (df['len\n(A)'] == len_i) &
        (df['sigma\n(A)'] == sigma_i) &
        (df['temp\n(k)'] == temp_i)
    ].index
    # 读取文件
    return df.at[row_number[0], item_i]
    # 项目名称对照表
    '''
        def __init__(self,):
            self.label_dict = {
                            "pore_radius\n(A)":                                                                                                            (r"Pore diameter ${D}_{\mathrm{GA}}$", r" $\mathrm{(nm)}$"), 
                            "pore_radius\n(A)(smoothed_by_Sigma)":                                                                                         (r"Pore diameter ${D}_{\mathrm{GA}}$", r" $\mathrm{(nm)}$"), 
                            "thickness_of_box\n(A)":                                                                                                       (r"GP thickness $t_{\mathrm{GA}}$", r" $\mathrm{ (\AA)}$"), 
                            "bond_density\n(unitless)":                                                                                                    (r"$C-C$ bond density ${\rho}_{C-C}$", r" (unitless)"),
                            "specific_surface_area\n(m^2/g)":                                                                                              (r"Specific surface area $S_{\mathrm{GA}}$", r" $\mathrm{(m^{2}/g)}$"), 
                            "specific_surface_area\n(m^2/g)(smoothed_by_Pore radius)":                                                                     (r"Specific surface area $S_{\mathrm{GA}}$", r" $\mathrm{(m^{2}/g)}$"), 
                            "surface_area\n(A^2)":                                                                                                         (r"Surface area $A_{\mathrm{GA}}$", r" $\mathrm{({\AA}^{2})}$"), 
                            "porosity\n(unitless)":                                                                                                        (r"Porosity $\phi_{\mathrm{GA}}$", r" ($\%$)"),
                            "porosity\n(unitless)(smoothed_by_Pore radius)":                                                                               (r"Porosity $\phi_{\mathrm{GA}}$", r" ($\%$)"),
                            "density\n(g/cm^3)":                                                                                                           (r"Density ${\rho}_{\mathrm{GA}}$",r" $\mathrm{ (mg/{cm}^{3})}$"),
                            'sigma\n(A)':                                                                                                                  (r"Sigma $\sigma$",r" $\mathrm{ (\AA)}$"),
                            'Diffusivity_EA_filtrates\n_in_membrane\n(m^2/s)':                                                                             (r"Diffusivity $\mathcal{D}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_TA_filtrates\n_in_membrane\n(m^2/s)':                                                                             (r"Diffusivity $\mathcal{D}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_TA_filtrates\n_in_membrane\n(m^2/s)(averaged)':                                                                   (r"Diffusivity $\mathcal{D}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_TA_filtrates\n_in_membrane\n(m^2/s)(averaged)(smoothed_by_Surface area)':                                         (r"Diffusivity $\mathcal{D}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_TA_filtrates\n_in_membrane\n(m^2/s)(averaged)(smoothed_by_Surface area)(smoothed_by_Pore radius)':                (r"Diffusivity $\mathcal{D}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_Knudsen_PoreRadius\n(m^2/s)':                                                                                     (r"Knudsen diffusivity $\mathcalD}_{Kn}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'Diffusivity_Knudsen_PoreRadius\n(m^2/s)Diffusivity_Knudsen_PoreRadius\n(m^2/s)(smoothed_by_Pore radius)':                     (r"Knudsen diffusivity $\mathcalD}_{Kn}$",r" $\mathrm{ (m^{2}/s)}$"),
                            'flux\n(L/m^2/h)':                                                                                                             (r"Flux $\mathcal{J}$",r" $\mathrm{(L/m^{2}/h)}$"),
                            'flux\n(L/m^2/h)(smoothed_by_Sigma)':                                                                                          (r"Flux $\mathcal{J}$",r" $\mathrm{(L/m^{2}/h)}$"),
                            'flux\n(L/m^2/h)(smoothed_by_Pore radius)':                                                                                    (r"Flux $\mathcal{J}$",r" $\mathrm{(L/m^{2}/h)}$"),
                            'flux\n(L/m^2/h)(smoothed_by_Pore radius)(smoothed_by_Tortuosity)':                                                            (r"Flux $\mathcal{J}$",r" $\mathrm{ (L/m^{2}/h)}$"),
                            'flux\n(L/m^2/h)(smoothed_by_Pore radius)(smoothed_by_Tortuosity)(smoothed_by_Porosity)':                                      (r"Flux $\mathcal{J}$",r" $\mathrm{ (L/m^{2}/h)}$"),
                            'tortuosity\n(unitless)':                                                                                                      (r"Tortuosity ${\tau}_{\mathrm{GA}}$", r" (unitless)"),
                            'tortuosity\n(unitless)(smoothed_by_Pore radius)':                                                                             (r"Tortuosity ${\tau}_{\mathrm{GA}}$", r" (unitless)"),
                            'tortuosity\n(unitless)(smoothed_by_Surface area)':                                                                            (r"Tortuosity ${\tau}_{\mathrm{GA}}$", r" (unitless)"),
                            'CN\n_cutoff':                                                                                                                 (r"$\#$ of $C$ atoms nearby $n_{C}$", r" (unitless)"),
                            'CN\n_1/2_peak':                                                                                                               (r"$\#$ of $C$ atoms nearby $n_{C}$", r" (unitless)"),
                            'CN\n_cutoff(smoothed_by_Pore radius)':                                                                                        (r"$\#$ of $C$ atoms nearby $n_{C}$", r" (unitless)"),
                            }
            self.color_lsit()
    '''



# read the data fule
class LmpGP_file_reader():
    def __init__(self,len,sigma,temp):
        ## terminology of parameters
        self.sigma = sigma
        self.len   = len
        self.temp  = temp
        ## working path 
        self.root           =  "../../Database"   
        self.As_root        =  f"{self.root}/OUTPUT/assemblys_pretreated"
        self.Di_root        =  f"{self.root}/OUTPUT/assemblys_distillated/len_{self.len}"  
        self.folder         =  f"len_{self.len}_sigma_{self.sigma}_{self.temp}"
        self.dump_file      =  f"dump_di"
        self.data_file      =  f"data.3_len_{self.len}_sigma_{self.sigma}_pretreated"
        self.temp_dump_file =  f"./merged_dump_file.txt"   
        ## auto excecuting
        #self.process_AGMD_dump_file()
        self.read_GP_data_file()
    def read_GP_data_file(self):
        ### initialize
        self.GP_atom_count, self.GP_box_bounds_info, self.GP_atom_data = None, {}, []
        file_path = os.path.join(f"{self.As_root}/{self.data_file}")
        atom_section = False

        ### extract info 
        with open(file_path, 'r') as file:
            for line in file:

                ## Get the number of atoms
                if line.strip().endswith("atoms"):
                    self.GP_atom_count = int(line.split()[0])

                ## Get the box info
                # boudnds
                if 'xlo xhi' in line:
                    self.GP_box_bounds_info['x'] = list(map(float, line.strip().split()[:2]))
                if 'ylo yhi' in line:
                    self.GP_box_bounds_info['y'] = list(map(float, line.strip().split()[:2]))
                if 'zlo zhi' in line:
                    self.GP_box_bounds_info['z'] = list(map(float, line.strip().split()[:2]))

                ## Read atom data
                # Check title
                if 'Atoms' in line and '#' in line:
                    atom_section = True
                    continue
                # read
                if atom_section and line.strip():
                    parts = line.split()
                    if len(parts) >= 6:  # valid line Ensured
                        atom_id, atom_type, x, y, z = int(parts[0]), int(parts[1]), float(parts[4]), float(parts[5]), float(parts[6])
                        self.GP_atom_data.append([atom_id, atom_type, x, y, z])
            
            # calculate side length
            x_length = float(self.GP_box_bounds_info['x'][1]) - float(self.GP_box_bounds_info['x'][0])
            y_length = float(self.GP_box_bounds_info['y'][1]) - float(self.GP_box_bounds_info['y'][0])
            z_length = float(self.GP_box_bounds_info['z'][1]) - float(self.GP_box_bounds_info['z'][0])
            self.GP_box_lengths_info = np.array([x_length, y_length, z_length])
            
        ### post-treatment
        ## convert
        self.GP_atom_data = np.array(self.GP_atom_data) 
        ## sorting
        # 首先按照 type 排序，然后在 type 相同的情况下按 ID 排序
        if hasattr(self, 'GP_atom_data') and self.GP_atom_data.size > 0:
            # 获取排序后的索引
            sorted_indices = np.lexsort((self.GP_atom_data[:, 0], self.GP_atom_data[:, 1]))  # 先按 ID，再按 type 排序
            # 应用排序
            self.GP_atom_data = self.GP_atom_data[sorted_indices]   
               
    def process_AGMD_dump_file(self):  
        ## Find dump file 
        path = f"{self.Di_root}/{self.folder}"
        # list all files
        file_list = os.listdir(path)
        # open folders
        file_list = [f for f in file_list if os.path.isfile(os.path.join(path, f))]
        # select dump file only
        file_list = [f for f in file_list if self.dump_file in f]
        # sort by time label
        file_list = [f.replace(f"{self.dump_file}_", "") for f in file_list]
        file_list.sort()
        file_list = [f"{path}/{self.dump_file}_{f}" for f in file_list]

        ## Combine dump files into one
        output_file_path = os.path.join("./merged_dump_file.txt")
        with open(output_file_path, 'w') as output_file:
            for file_name in file_list:
                with open(os.path.join(file_name), 'r') as input_file:
                    output_file.write(input_file.read())
        
        ### Collect basic info 
        ## Comput TIMESTEP positions list
        with open(self.temp_dump_file, 'r') as file:
            for line in file:
                ## box info
                if 'ITEM: BOX BOUNDS' in line:
                    self.box_bounds_info = {}
                    self.box_bounds_info["x"] = list(map(float, next(file).strip().split()))
                    self.box_bounds_info["y"] = list(map(float, next(file).strip().split()))
                    self.box_bounds_info["z"] = list(map(float, next(file).strip().split()))
                    x_length = float(self.box_bounds_info['x'][1]) - float(self.box_bounds_info['x'][0])
                    y_length = float(self.box_bounds_info['y'][1]) - float(self.box_bounds_info['y'][0])
                    z_length = float(self.box_bounds_info['z'][1]) - float(self.box_bounds_info['z'][0])
                    self.box_lengths_info = np.array([x_length, y_length, z_length])
                    break
        ## build timesteps list
        self.timestep_positions, self.timestep_list = {}, []
        timestep_checker = 200000
        with open(self.temp_dump_file, 'r') as file:
            line_index = -1
            while True:
                line = file.readline()
                line_index += 1
                if not line:
                    break  # 到达文件末尾
                if 'ITEM: TIMESTEP' in line:
                    # 读取时间步的值
                    timestep = int(file.readline().strip())
                    line_index += 1
                    if timestep != timestep_checker:
                        # 记录当前时间步的起始行号
                        self.timestep_positions[timestep] = line_index 
                        # 将步数记入列表 
                        self.timestep_list.append(timestep) 
                        timestep_checker = timestep        
    def read_single_timestep_for__AGMD_dump_file(self, timestep):
        #### 打开文件并开始读取
        file_path = os.path.join(self.temp_dump_file)
        with open(file_path, 'r') as file:
            ### Jump to timestep
            file.seek(0)  # 开始时先回到文件开头
            for _ in range(self.timestep_positions[timestep]):  # 跳过line_number - 1行
                file.readline()
            ## Current number of atoms
            while True:
                line = next(file)
                if 'ITEM: NUMBER OF ATOMS' in line:
                    self.atom_count = int(next(file).strip())
                    break
            ## Skip
            while True:
                line = next(file)
                if 'ITEM: ATOMS' in line:
                    break

            ### Extract atoms info
            atom_data = []
            for _ in range(self.atom_count):
                # 读取并分割行
                atom_info_str = next(file).strip().split()
                # 分别处理不同类型的数据
                # ITEM: ATOMS id type xs ys zs
                atom_info = [int(atom_info_str[0]), int(atom_info_str[1])] + [float(x) for x in atom_info_str[2:]]
                atom_data.append(atom_info)

        ### Post treatment
        # convert to NumPy matrix
        self.atom_data = np.array(atom_data) 
        ## de-normalize
        if hasattr(self, 'box_lengths_info') and self.atom_data.size > 0:
            for i in range(self.atom_data.shape[0]):
                # xs 
                self.atom_data[i, 2] = self.atom_data[i, 2] * self.box_lengths_info[0] + float(self.box_bounds_info["x"][0])
                # ys 
                self.atom_data[i, 3] = self.atom_data[i, 3] * self.box_lengths_info[1] + float(self.box_bounds_info["y"][0])
                # zs
                self.atom_data[i, 4] = self.atom_data[i, 4] * self.box_lengths_info[2] + float(self.box_bounds_info["z"][0])
        ## sorting
        # 首先按照 type 排序，然后在 type 相同的情况下按 ID 排序
        if hasattr(self, 'atom_data') and self.atom_data.size > 0:
            # 获取排序后的索引
            sorted_indices = np.lexsort((self.atom_data[:, 0], self.atom_data[:, 1]))  # 先按 ID，再按 type 排序
            # 应用排序
            self.atom_data = self.atom_data[sorted_indices]
    def build_dump_atoms_info_list__for__AGMD_dump_file(self):
        # 初始化一个列表来存储每个时间步的原子信息矩阵
        self.atom_data_list = []
        # 遍历每个时间步
        for timestep in self.timestep_positions.keys():
            # 读取该时间步的原子信息
            self.read_single_timestep_for__AGMD_dump_file(timestep)
            # 将原子信息矩阵添加到列表中
            self.atom_data_list.append(self.atom_data)
            print(timestep)
    def check_atom_data_list(self):
        self.timestep_list
        self.box_bounds_info
    def check_atom_data_list(self, output_filename="check_atom_data_list.txt"):
        """
        将atom_data_list的内容输出为LAMMPS dump文件格式
        :param output_filename: 输出文件名
        """
        with open(output_filename, 'w') as dump_file:
            for i, atom_data in enumerate(self.atom_data_list):
                # 写入时间步
                dump_file.write("ITEM: TIMESTEP\n")
                dump_file.write(f"{self.timestep_list[i]}\n")
                # 写入原子数量
                dump_file.write("ITEM: NUMBER OF ATOMS\n")
                dump_file.write(f"{len(atom_data)}\n")
                # 写入盒子边界
                dump_file.write("ITEM: BOX BOUNDS pp pp ff\n")
                for bound in ['x', 'y', 'z']:
                    dump_file.write(f"{self.box_bounds_info[bound][0]} {self.box_bounds_info[bound][1]}\n")
                # 写入原子信息
                dump_file.write("ITEM: ATOMS id type x y z\n")
                for atom in atom_data:
                    dump_file.write(" ".join(map(str, atom)) + "\n")

