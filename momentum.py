import uproot
import numpy as np
import os
import vector
eps = 1e-2
def is_pt_same(data):
    file = uproot.open(data)
    tree = file['tree']

    jets = {}

    #v vsakem jet-u so delci in njihove gibalne količine
    components = ['px', 'py', 'pz']
    
    for component in components:
        part = 'part_'+component
        branch = tree[part]
        # jets['jets_'+component] = [jet_px[part] for jet_px in branch.arrays().tolist()] #tale dela ziher
        jets[component] = [jet_px for jet_px in branch.array().tolist()]
    
    branch = tree['jet_pt']
    jet_pts1 = sorted(branch.array())

    branch = tree['jet_eta']
    jet_eta1 = sorted(branch.array())

    jet_eta2 = []
    jet_pts2 = []
    
    for i in range(len(jet_pts1)):
        ###transversalna gibalna kolicina
        jet_pt = np.hypot(sum(jets['px'][i]), sum(jets['py'][i]))
        jet_pts2.append(jet_pt)
        ###rapidnost
        pz = sum(jets['pz'][i])
        p = np.hypot(jet_pt,pz)
        # eta=np.arctanh(pz/np.hypot(jet_pt, pz))
        jet_eta2.append(0.5*np.log((p+pz)/(p-pz)))

    jet_pts2.sort()
    jet_eta2.sort()

    # print(*jet_eta2[:10])
    # print(*jet_eta1[:10])
    # exit()
    for i in range(len(jet_pts1)):
        if abs(jet_pts1[i] - jet_pts2[i]) > eps:
            file.close()
            print(i, jet_pts1[i] - jet_pts2[i], jet_pts1[i], jet_pts2[i])
            return False
        if abs(jet_eta1[i] - jet_eta2[i]) > eps:
            file.close()
            print(i, jet_eta1[i] - jet_eta2[i], jet_eta1[i], jet_eta2[i])
            return False
    
    file.close()
    return True
    #isci za vse delce v jetu skupno gibalno, nakonc bi moglo bit isto 

prefix = "./data/"
files = os.listdir(prefix)
# files = ["TTBar_100.root"]
for i in files:
    print("{} is {}".format(i ,"OK" if is_pt_same(prefix+i) else "FAIL"))
    pass