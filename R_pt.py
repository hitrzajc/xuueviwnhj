import uproot
import numpy as np
import os
import matplotlib.pyplot as plt
eps = 1e-2
def geteta(p,pz):
    return 0.5*np.log((p+pz)/(p-pz))
def getphi(py,px):
    return np.arctan2(py,px)

def moment(p,hist):
    ans = 0
    for i in range(len(hist[0])):
        ans += pow(hist[1][i],p)*hist[0][i]
    return ans

def cmoment(p,hist): #centralni moment
    ans = 0
    xbar = 0
    for i in range(len(hist[0])):
        xbar += hist[0][i]*hist[1][i]
    for i in range(len(hist[0])):
        ans += pow(hist[1][i]-xbar,p)*hist[0][i]
    return ans
    
    
def risi(data):
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
    
    R = [[] for i in range(5)]
    Pt = [[] for i in range(5)]
    njets = [0]*5
    for k in range(len(jets['px'])):
        px = sum(jets['px'][k])
        py = sum(jets['py'][k])
        pz = sum(jets['pz'][k])
        pt = 0 #tale bo za normalizacijo
        for i in range(len(jets['px'][k])):
            pxi = jets['px'][k][i]
            pyi = jets['py'][k][i]
            pt += np.hypot(pxi,pyi)
        ptt = np.hypot(px,py)
        p = np.hypot(ptt,pz)
        phi = getphi(py,px)
        eta = geteta(p,pz)
        index = int(ptt)//100-5
        tmpR = R[index]
        tmpPt = Pt[index]
        njets[index] += 1
        for i in range(len(jets['px'][k])):
            pxi = jets['px'][k][i]
            pyi = jets['py'][k][i]
            pzi = jets['pz'][k][i]
            pti = np.hypot(pxi,pyi)
            pi = np.hypot(pti,pzi)
            etai = geteta(pi,pzi)
            phii = getphi(pyi,pxi)
            dphi = phii - phi
            if(dphi < -np.pi): dphi+=2*np.pi
            if(dphi > +np.pi): dphi-=2*np.pi
            
            deta = etai - eta
            dRi = np.sqrt(dphi**2 + deta**2)
            if dRi>1:
                continue
            tmpR.append(dRi)
            # tmpPt.append((pt,pti))
            tmpPt.append((pt,pi))
    for i in range(5):
        a = plt.hist(R[i], bins=100,weights=[pti/pt/njets[i] for pt,pti in Pt[i]], alpha=0.7, color='blue')
        textstr=""
        for j in range(1,5):
            textstr+=r'$M_{%i}={%f}$' % (j,moment(j,a)) + '\n'
        
        plt.text(0.95, 0.95, textstr, 
         transform=plt.gca().transAxes, 
         verticalalignment='top', 
         horizontalalignment='right', 
         fontsize=12, 
         bbox=dict(facecolor='white', alpha=0.5))

        # print(sum(a[0]))
        plt.title('Histogram for {}00 to {}00 GeV'.format(i+5,i+6))
        plt.xlabel('$R$')
        plt.ylabel('$P - probability$')
        plt.xlim(0,1)

        plt.savefig("{}.pdf".format(i+1))
        plt.clf()
        


prefix = "./data/"
files = os.listdir(prefix)
files = ["TTBar_100.root"]
for i in files:
    # print("{} is {}".format(i ,"OK" if is_pt_same(prefix+i) else "FAIL"))
    risi(prefix+i)
    pass