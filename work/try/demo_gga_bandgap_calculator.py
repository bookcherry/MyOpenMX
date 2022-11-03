from pymatgen.ext.matproj import MPRester
import os
import pandas as pd
import subprocess
import numpy as np
import csv

API_KEY = 'hogehoge'
#hogehogeの部分にMaterials ProjectのAPI_Keyをいれる。

mat_list_name = "mat_list.txt"
xc_type = "GGA-PBE"
unit = 27.2113845


def check_string(search_str):
  with open(f"demo_GGA-PBE_result.csv") as temp_f:
    datafile=temp_f.readlines()
    for temp_line in datafile:
      if search_str in temp_line:
        return True
    return False

os.makedirs(f'./demo_GGA-PBE', exist_ok = True)

os.chdir(f"./demo_GGA-PBE")

with open("../demo_mat_list.csv", "r") as f:
  reader = csv.reader(f)
  for row in reader:
    line = row[0]
    experimental_band_gap = float(row[1])
    benchmark_gga_gap = float(row[2])
    if len(line)!=0:
      with MPRester(API_KEY) as m:
          # results = ,.query(properties =)
          data = m.query(criteria={"task_id":line},properties=["cif", "band_gap", "pretty_formula"])
        
          mat_name = data[0]["pretty_formula"]
          if check_string(f"{mat_name}_{line}"):
            continue        
          os.makedirs(f"./{mat_name}_{line}", exist_ok = True)
          with open(f"./{mat_name}_{line}/{mat_name}_{line}.cif", "w") as g:
              try:
                  g.write(data[0]["cif"])
              except IndexError:
                  print("Error")
          ref_dft_band_gap = data[0]["band_gap"]
          
  
          # os.chdir("/home/k-tsukamoto/Tool/openmx3.9/work/try")
          os.chdir(f"./{mat_name}_{line}")
    
          
          subprocess.run(['python3', '/home/ozaki-lab/cif2input/cif2input.py', f"{mat_name}_{line}.cif"])
          
          
          input_file = "openmx.dat"
          
          with open(input_file, encoding="cp932") as f:
              data_lines = f.read()
          
          # 文字列置換
          data_lines = data_lines.replace("System.Name          openmx", f"System.Name          {mat_name}_{line}")
          data_lines = data_lines.replace("scf.XcType               GGA-PBE", f"scf.XcType               {xc_type}")
          data_lines = data_lines.replace("data.path     ../DFT_DATA19/", "data.path     ../../../../DFT_DATA19")
          data_lines = data_lines.replace("Band.dispersion              off", "Band.dispersion              on")
          
          
          # 同じファイル名で保存
          with open(input_file, mode="w", encoding="cp932") as f:
              f.write(data_lines)
          
          subprocess.run(['mpirun', '-np', '4', '../../../openmx', input_file, '|', 'tee', f"{mat_name}_{line}.std"])
          
          command = f"grep 'g1 =' '{mat_name}_{line}.out'" + "| awk '{print $3}'"
          g1 = float(subprocess.run(command, shell = True, capture_output= True, text = True).stdout)
          command = f"grep 'g2 =' '{mat_name}_{line}.out'" + "| awk '{print $3}'"
          g2 = float(subprocess.run(command, shell = True, capture_output= True, text = True).stdout)
          command = f"grep 'KSGap =' '{mat_name}_{line}.out'" + "| awk '{print $3}'"
          KSGapFromKMesh = float(subprocess.run(command, shell = True, capture_output= True, text = True).stdout)
          command = f"grep 'Di_Pe =' '{mat_name}_{line}.out'" + "| awk '{print $3}'"
          Di_Pe = float(subprocess.run(command, shell = True, capture_output= True, text = True).stdout)
          
          subprocess.run(['../../../bandgnu13', f'{mat_name}_{line}.Band'])
            
          band_path = f'{mat_name}_{line}.BANDDAT1'
            
          count = 0
            
          command = "awk  '/<Band.kpath/,/Band.kpath>/' " + f"{input_file}" + " | wc -l"        
         
            
          num_spaces = subprocess.run(command, shell = True, capture_output= True, text = True).stdout
 
          num_spaces = int(num_spaces)-2
            
          command = "awk 'BEGIN{a=0} /<Atoms.SpeciesAndCoordinates/,/Atoms.SpeciesAndCoordinates>/{a+=$6+$7} END{print a}' " + f"{input_file}"        
            
          num_valence_electrons = subprocess.run(command, shell = True, capture_output= True, text = True).stdout
            
          num_orbitals = int(num_valence_electrons)/2
            
          first_breakline = num_spaces-1
          
          homo_read_line = 2*(num_spaces*(num_orbitals-1)+first_breakline)
          lumo_read_line = 2*(num_spaces*(num_orbitals)+first_breakline)
          secondlumo_read_line = 2*(num_spaces*(num_orbitals+1)+first_breakline)
          homo_rows = []
          lumo_rows = []
          with open(band_path) as f:
              reader = csv.reader(f)
              for row in reader:
                  if row == []:
                      count += 1
                  if row != []:
                      if homo_read_line <= count < lumo_read_line:
                          row = row[0].split()
                          homo_rows.append(row)
                      if lumo_read_line <= count < secondlumo_read_line:
                          row = row[0].split()
                          lumo_rows.append(row)
                        
          chemP = float(subprocess.run("grep 'Chemical Potential' " + f"*_{line}.out " + "| awk '{print $5}'", shell = True, capture_output= True, text = True).stdout)*unit
        #   print(chemP)
          
          
          homo_rows  = [[float(x) for x in y] for y in homo_rows]
          lumo_rows  = [[float(x) for x in y] for y in lumo_rows]
          
          homo_rows = np.array(homo_rows)
          
          homo_rows = homo_rows[:,1]
          
          lumo_rows = np.array(lumo_rows)
          
          lumo_rows = lumo_rows[:,1]
          
          lumo_energy = np.min(lumo_rows) + chemP
          homo_energy = np.max(homo_rows) + chemP
          
          KSGapFromBandStructure = lumo_energy - homo_energy
          
          os.chdir("..")
          
          with open(f"demo_GGA-PBE_result.csv", "a") as f:
            f.write(f"{mat_name}_{line},{g1},{g2},{Di_Pe},{KSGapFromKMesh},{KSGapFromBandStructure},{benchmark_gga_gap},{experimental_band_gap}\n")
           

   
      
      
      
      
      
# execute cif2input

# modify the openmx.dat

# execute the calculation

# calculate the band gap

# added the band gap file (local and global) (experimental and calculation( with xc functional))

  
  
  