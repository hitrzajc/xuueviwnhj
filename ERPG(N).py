import uproot
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def geteta(p,pz):
    return 0.5*np.log((p+pz)/(p-pz))
def getphi(py,px):
    return np.arctan2(py,px)

def moment(p,R,f):
    ans = 0
    for i in range(len(R)):
        ans += pow(R[i],p)*f[i]
    return ans

def cmoment(p,R,f): #centralni moment
    ans = 0
    xbar = 0
    for i in range(len(R)):
        xbar += R[i]*f[i]
    for i in range(len(R)):
        ans += pow(R[i]-xbar,p)*f[i]
    return ans

data = "./data/TTBar_100.root"
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

data = []
particles = [[] for i in range(len(jets['px']))]
def calc_theta(p1,p2):
    phi1 = p1[2]
    phi2 = p2[2]
    eta1 = p1[1]
    eta2 = p2[1]
    dphi = phi1 - phi2
    if(dphi < -np.pi): dphi+=2*np.pi
    if(dphi > +np.pi): dphi-=2*np.pi
    return np.sqrt(dphi**2 + (eta1-eta2)**2)

def f(k, i, n, pti, theta, pt):
    if n==0:
        if theta>1:
            return 0
        data.append((theta,2*pti/pow(pt,N)))
        return 0
    # out_pt = 0
    # out_theta = 0
    if n==N:
        for j in range(0,len(particles[k])):
            f(k,j,n-1,
              pti*particles[k][j][0], 
              0,pt)
    else:
        for j in range(i+1,len(particles[k])):
            f(k,j,n-1,
              pti*particles[k][j][0], 
              calc_theta(particles[k][i],particles[k][j]),pt)
            # out_pt += tmp_pt
            # out_theta += tmp_theta
        # return (out_pt, out_theta)
N=2
for k in range(len(jets['px'])):
    if k>5000:break
    print(k,end='\r')
    px = sum(jets['px'][k])
    py = sum(jets['py'][k])
    pz = sum(jets['pz'][k])
    pt = 0 #tale bo za normalizacijo
    for i in range(len(jets['px'][k])):
        pxi = jets['px'][k][i]
        pyi = jets['py'][k][i]
        pt += np.hypot(pxi,pyi)
    # ptt = np.hypot(px,py)
    # p = np.hypot(ptt,pz)
    # phi = getphi(py,px)
    # eta = geteta(p,pz)
    for i in range(len(jets['px'][k])):
        pxi = jets['px'][k][i]
        pyi = jets['py'][k][i]
        pzi = jets['pz'][k][i]
        pti = np.hypot(pxi,pyi)
        pi = np.hypot(pti,pzi)
        etai = geteta(pi,pzi)
        phii = getphi(pyi,pxi)

        particles[k].append((pti,etai,phii))
    f(k,0,N,1,0,pt)
    # ptt/=pt*pt
    # data.append((2*ptt,theta))
    # N=2
    # for i in range(len(particles)):
    #     pti,etai,phii =  particles[i]
    #     for j in range(i+1,len(particles)):
    #         ptj,etaj,phij = particles[j]
    #         deta = etai - etaj
    #         dphi = phii - phij
    #         if(dphi < -np.pi): dphi+=2*np.pi
    #         if(dphi > +np.pi): dphi-=2*np.pi
    #         dtheta = np.sqrt(dphi**2 + deta**2)
    #         if(dtheta>1):continue
    #         data.append((dtheta,2*pti*ptj/(pt*pt)))

x = [i[0] for i in data]
weights = [i[1]/5000 for i in data]
# for i in weights:
plt.hist(x,bins=100, weights=weights)
plt.show()

