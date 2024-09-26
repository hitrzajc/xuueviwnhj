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
for k in range(len(jets['px'])):
    print("{}/{}".format(k,len(jets['px'])),end='\r')
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
    tmpR = []
    f = []
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
        # dRi = dphi**2 + deta**2
        
        if dRi>1:
            continue
        tmpR.append(dRi)
        f.append(pti/pt)
    
    data.append((pt,
                 moment(1,tmpR,f),
                 moment(2,tmpR,f),
                 moment(3,tmpR,f),
                 moment(4,tmpR,f)))
data.sort()
mind = min(data)[0]
# maxd = max(data)[0]
maxd = 1000
bsize = (maxd-mind)/100

c = mind
i=0
baskets = [[] for i in range(100)]
for d in data:
    p,m1,m2,m3,m4 = d
    if p>maxd: continue
    if c+(i+1)*bsize<p:
        i+=1
        c+=bsize
    baskets[i].append(d)

def get_avg(basket,i):
    ans=0
    for m in basket:
        ans+=m[i]
    return ans/len(basket)

# def get_error(basket,i):
#     sum = 0
#     ans = 0
#     for k in basket:
#         sum += k[i]
#     avg = sum/len(basket)
#     for k in basket:
#         ans += pow((k[i]-avg),2)
#     ans/=sum
#     return np.sqrt(ans/len(basket))
    # return None
def get_error(basket,i):
    sum = 0
    ans = 0
    for k in basket:
        sum += k[i]
    avg = sum/len(basket)
    for k in basket:
        ans += pow((k[i]-avg),2)
    ans = np.sqrt(ans/len(basket))
    return ans/np.sqrt(len(basket))*1.96
x = []
for basket in baskets:
    if len(basket) == 0:
        continue
    x.append(get_avg(basket,0))

def model(x,alpha,A):
    p0 = 600
    x = np.array(x)
    return A*np.float_power(p0/x,alpha)

for i in range(1,5):
    y=[]
    error = [] #s

    for basket in baskets:
        if len(basket) == 0:
            continue
        avg = get_avg(basket,i)
        y.append(avg)
        error.append(get_error(basket,i))
    # print(error)
    plt.errorbar(x, y, yerr=error, fmt='o', capsize=5, label='$M$', zorder=1)
    # plt.title('$M_{}(p_t)$ za centralne momente distribucije $R$'.format(i))
    plt.title('$M_{}(p_t)$ za momente distribucije $R$'.format(i))

    plt.xlabel('$p_t[GeV]$')
    plt.ylabel('$M_{}$'.format(i))

    filterx = []
    filtery = []
    filtere = []
    for k in range(len(x)):
        if(550<=x[k] and x[k]<=900):
            filterx.append(x[k])
            filtery.append(y[k])
            filtere.append(error[k])
    popt, pcov = curve_fit(model,filterx , filtery,sigma=filtere, absolute_sigma=True)
    alpha, A = popt
    alpha_err, A_err = np.sqrt(np.diag(pcov))  
    plt.plot(filterx, model(filterx, *popt), color='red',zorder =2,
              label=r'$f(p_t) = A\left(\frac{p_0}{p_t}\right)^\alpha$;$p_0=600[GeV]$, $A={%.3f}$, $\alpha={%.3f}\pm{%.3f}$' % (A,alpha,alpha_err))
    plt.legend()
    plt.savefig("CCM_{}(P_t).pdf".format(i))

    plt.clf()


