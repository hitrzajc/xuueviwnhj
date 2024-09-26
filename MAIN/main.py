import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import functions
from tqdm import tqdm

import concurrent.futures

# Specify the directory
directory = "/home/tadej/work/test_20M/"

# List files in the directory
files = os.listdir(directory)

names = {}

for file in files:
    name = file[:-9]
    num = file[-8:-5]
    if name not in names:
        names[name] = []
    names[name].append(num)


#MULTI THREADING
pool = concurrent.futures.ThreadPoolExecutor(max_workers=32)

##### ERPG
print("""
▓█████ ██▀███  ██▓███   ▄████ 
▓█   ▀▓██ ▒ ██▓██░  ██▒██▒ ▀█▒
▒███  ▓██ ░▄█ ▓██░ ██▓▒██░▄▄▄░
▒▓█  ▄▒██▀▀█▄ ▒██▄█▓▒ ░▓█  ██▓
░▒████░██▓ ▒██▒██▒ ░  ░▒▓███▀▒
░░ ▒░ ░ ▒▓ ░▒▓▒▓▒░ ░  ░░▒   ▒ 
 ░ ░  ░ ░▒ ░ ▒░▒ ░      ░   ░ 
   ░    ░░   ░░░      ░ ░   ░ 
   ░  ░  ░                  ░ 
                              
""")

bins = 100
K = 0 #num of jets in each file if 0 all jets
for name in names:
    # name = "HToGG"
    def F1(names,name,bins,K):
        nums = names[name]
        weight = np.zeros(bins)
        for num in tqdm(nums, desc=f"Processing Files for {name}", leave=True):
            file = directory + name + "_" + num + ".root"
            # print("WORK: {}".format(file))
            weight += functions.ERPG_2(file, K, bins)
        weight/=len(nums)
        plt.hist([i/bins for i in range(bins)],bins=bins, weights=weight)
        plt.xlabel("$\\vartheta$")
        plt.ylabel("$\\frac{p_ip_j}{p_t^2}$")
        plt.savefig("pdfs/" + "ERPG_2_" +name + ".pdf")
        plt.clf()
    pool.submit(F1, names, name, bins, K)
###########


###### M(P_T)
print("""
                                                          
@@@@@@@@@@      @@@  @@@@@@@             @@@@@@@  @@@     
@@@@@@@@@@@    @@@   @@@@@@@@            @@@@@@@   @@@    
@@! @@! @@!   @@!    @@!  @@@              @@!      @@!   
!@! !@! !@!  !@!     !@!  @!@              !@!       !@!  
@!! !!@ @!@  !!@     @!@@!@!               @!!       !!@  
!@!   ! !@!  !!!     !!@!!!                !!!       !!!  
!!:     !!:  !!:     !!:                   !!:       !!:  
:!:     :!:   :!:    :!:                   :!:      :!:   
:::     ::      ::    ::  :::::::::::::     ::     ::     
 :      :         :   :   :::::::::::::     :     :       
                                                          
""")
for M in range(1,5):
    for name in names:break #resets variable "name"
    for name in tqdm(names, desc=f"Processing {M}/4 {name}", leave=True):
        def F2(names,name,M):
            nums = names[name]
            file = directory + name + "_" + nums[0] + ".root"
            functions.M_R(name,file,M=M)
        pool.submit(F2,names,name,M)


######## p_t(R)
print("""
██████╗ ████████╗ ██╗██████╗ ██╗ 
██╔══██╗╚══██╔══╝██╔╝██╔══██╗╚██╗
██████╔╝   ██║   ██║ ██████╔╝ ██║
██╔═══╝    ██║   ██║ ██╔══██╗ ██║
██║███████╗██║   ╚██╗██║  ██║██╔╝
╚═╝╚══════╝╚═╝    ╚═╝╚═╝  ╚═╝╚═╝ 
            
""")
bins = 100
N = 4
bin_siz = 1/bins
for name in names:
    def F3(bins,N,names,name):
        functions.pt_R_plot(names,name,directory,bins,"pt_R",N)
    def F4(bins,N,names,name):
        functions.pt_R_plot(names,name,directory,bins,"pt_R2",N)
    pool.submit(F3,bins,N,names,name)
    pool.submit(F4,bins,N,names,name)

    
pool.shutdown(wait=True)
