import uproot
import numpy as np
import matplotlib.pyplot as plt


def cmoment(p,hist): #centralni moment
    ans = 0
    xbar = 0
    for i in range(len(hist[0])):
        xbar += hist[0][i]*hist[1][i]
    for i in range(len(hist[0])):
        ans += pow(hist[1][i]-xbar,p)*hist[0][i]
    return ans

data = "./data/TTBar_100.root"
# myFile = ROOT.TFile.Open(data, "RECREATE")
# print(myFile.px())

file = uproot.open(data)
# print(file['tree'].values())
# print(file['tree'].axis().edges())
tree = file['tree']
eps = 1e-2
def geteta(p,pz):
    return 0.5*np.log((p+pz)/(p-pz))
def getphi(py,px):
    return np.arctan2(py,px)

components = ['px','py','pz']
jets = {}
for component in components:
    part = 'part_'+component
    branch = tree[part]
    # jets[component] = [jet_px[part] for jet_px in branch.arrays().tolist()] #tale dela ziher
    jets[component] = [jet for jet in branch.array().tolist()]
    
R = [[] for i in range(5)]
Pt = [[] for i in range(5)]
for k in range(5,6):
    px = sum(jets['px'][k])
    py = sum(jets['py'][k])
    pz = sum(jets['pz'][k])
    pt = 0
    for i in range(len(jets['px'][k])):
        pxi = jets['px'][k][i]
        pyi = jets['py'][k][i]
        pt += np.hypot(pxi,pyi)
    p = np.hypot(pt,pz)
    phi = getphi(py,px)
    eta = geteta(p,pz)
    tmpR = R[int(pt)//100-5]
    tmpPt = Pt[int(pt)//100-5]
    BBB = int(pt)//100-5
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
        tmpPt.append((pt,pti))

i=BBB
# weights=
# print(*R[i])
# print(*Pt[i])
s = 0
for pt,pti in Pt[i]:
    s += (pti/pt)
a = plt.hist(R[i], bins=10,weights=[pti/pt for pt,pti in Pt[i]], alpha=0.7, color='blue')
textstr=""
for j in range(1,5):
    textstr+=r'$M_{%i}={%f}$' % (j,cmoment(j,a)) + '\n'

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
plt.show()
    
