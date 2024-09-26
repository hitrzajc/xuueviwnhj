import numpy as np
import uproot
from tqdm import tqdm
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

def ERPG_2(file, K, bins, N=2):
    file = uproot.open(file)
    tree = file['tree']
    jets = {}
    #v vsakem jet-u so delci in njihove gibalne količine
    components = ['px', 'py', 'pz']
    if K == 0:
        K = len(tree['part_px'].array())
    for component in components:
        part = 'part_'+component
        branch = tree[part]
  
        jets[component] = [jet_px for jet_px in branch.array().tolist()]


    weights = np.zeros(bins)
    # def calc_theta(p1,p2):
    #     phi1 = p1[2]
    #     phi2 = p2[2]
    #     eta1 = p1[1]
    #     eta2 = p2[1]
    #     dphi = phi1 - phi2
    #     if(dphi < -np.pi): dphi+=2*np.pi
    #     if(dphi > +np.pi): dphi-=2*np.pi
    #     return np.sqrt(dphi**2 + (eta1-eta2)**2)

    bin_siz = 1/bins
    def query(theta, pt):
        bin_idx = int(theta/bin_siz)
        weights[bin_idx] += pt
    
    if K == 0:
        tr = tqdm(range(len(jets['px'])),leave=False)
    else:
        tr = tqdm(range(len(jets['px'])),total=K,leave=False)
        
    for k in tr:
        if k>K:break
        # print("{}/{}".format(k,K),end='\r')
        pt = 0 #tale bo za normalizacijo
        for i in range(len(jets['px'][k])):
            pxi = jets['px'][k][i]
            pyi = jets['py'][k][i]
            pt += np.hypot(pxi,pyi)
        
        particles = []
        for i in range(len(jets['px'][k])):
            pxi = jets['px'][k][i]
            pyi = jets['py'][k][i]
            pzi = jets['pz'][k][i]
            pti = np.hypot(pxi,pyi)
            pi = np.hypot(pti,pzi)
            etai = geteta(pi,pzi)
            phii = getphi(pyi,pxi)

            particles.append((pti,etai,phii))

        for i in range(len(particles)):
            pti,etai,phii =  particles[i]
            for j in range(i+1,len(particles)):
                # if j==i: continue
                ptj,etaj,phij = particles[j]
                deta = etai - etaj
                dphi = phii - phij
                if(dphi < -np.pi): dphi+=2*np.pi
                if(dphi > +np.pi): dphi-=2*np.pi
                dtheta = np.sqrt(dphi**2 + deta**2)
                if(dtheta>1):continue
                query(dtheta,2*pti*ptj/(pt*pt))

    weights/=K
    return weights

def M_R(name,file, bins=100, M=1):
    tree = uproot.open(file)['tree']
    jets = {}
    #v vsakem jet-u so delci in njihove gibalne količine
    components = ['px', 'py', 'pz']

    for component in components:
        part = 'part_'+component
        branch = tree[part]
        # jets['jets_'+component] = [jet_px[part] for jet_px in branch.arrays().tolist()] #tale dela ziher
        jets[component] = [jet for jet in branch.array().tolist()]

    data=[]
    for k in tqdm(range(len(jets['px'])),desc=f"{file}",leave=False):
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
            # dRi = np.sqrt(dphi**2 + deta**2)
            dRi = dphi**2 + deta**2

            if dRi>1:continue
            tmpR.append(dRi)
            f.append(pti/pt)

        data.append((pt,moment(M,tmpR,f)))
    data.sort()
    mind = min(data)[0]
    # maxd = max(data)[0]
    maxd = 1000
    bsize = (maxd-mind)/100

    c = mind
    i=0
    baskets = [[] for i in range(100)]
    for d in data:
        p,m = d
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

    def get_error(basket,i):
        sum = 0
        ans = 0
        for k in basket:
            sum += k[i]
        avg = sum/len(basket)
        for k in basket:
            ans += pow((k[i]-avg),2)
        ans/=sum
        return np.sqrt(ans/len(basket))
        # return None

    x = []
    for basket in baskets:
        if len(basket) == 0:
            continue
        x.append(get_avg(basket,0))

    def model(x,alpha,A):
        p0 = 600
        x = np.array(x)
        return A*np.float_power(p0/x,alpha)
    
    y=[]
    error = [] #s
    for basket in baskets:
        if len(basket) == 0:
            continue
        avg = get_avg(basket,1)
        y.append(avg)
        error.append(get_error(basket,1))
    plt.errorbar(x, y, yerr=error, fmt='o', capsize=5, label='$M$', zorder=1)
    # plt.title('$M_{}(p_t)$ za centralne momente distribucije $R$'.format(i))
    plt.title('$M_{}(p_t)$ za momente distribucije $R^2$'.format(M))
    plt.xlabel('$p_t[GeV]$')
    plt.ylabel('$M_{}$'.format(M))
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
    plt.savefig("pdfs/M_{}(P_t)_".format(M)+name+".pdf")
    plt.clf()

def pt_R(file , bins): #za pt(R)
    tree = uproot.open(file)['tree']
    jets = {}

    #v vsakem jet-u so delci in njihove gibalne količine
    components = ['px', 'py', 'pz']
    
    for component in components:
        part = 'part_'+component
        branch = tree[part]
        # jets['jets_'+component] = [jet_px[part] for jet_px in branch.arrays().tolist()] #tale dela ziher
        jets[component] = [jet_px for jet_px in branch.array().tolist()]
    
    weights = np.zeros(bins)
    # pt_mx = 1000
    # pt_mn = 500

    bin_size = 1/bins
    for k in tqdm(range(len(jets['px']))):
        # if k>5000:break #tole zbirsi za full scale
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
            # dRi = dphi + deta
            if dRi>1:
                continue
            index = int(dRi/bin_size)
            weights[index] += pti/pt
    #tole zamenjaj za ful scale
    weights /= len(jets['px'])
    # weights/=5000
    return np.array(weights)

def pt_R2(file , bins): #za pt(R^2)
    tree = uproot.open(file)['tree']
    jets = {}

    #v vsakem jet-u so delci in njihove gibalne količine
    components = ['px', 'py', 'pz']
    
    for component in components:
        part = 'part_'+component
        branch = tree[part]
        # jets['jets_'+component] = [jet_px[part] for jet_px in branch.arrays().tolist()] #tale dela ziher
        jets[component] = [jet_px for jet_px in branch.array().tolist()]
    
    weights = np.zeros(bins)
    # pt_mx = 1000
    # pt_mn = 500

    bin_size = 1/bins
    for k in tqdm(range(len(jets['px']))):
        # if k>5000:break #tole zbirsi za full scale
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
            # dRi = np.sqrt(dphi**2 + deta**2)
            dRi = dphi**2 + deta**2
            if dRi>1:
                continue
            index = int(dRi/bin_size)
            weights[index] += pti/pt
    #tole zamenjaj za ful scale
    weights /= len(jets['px'])
    # weights/=5000
    return weights

def pt_R_plot(names, name, directory, bins, function_name, N=4):
    bin_siz = 1/bins
    x_data = [i*bin_siz for i in range(bins)]
    nums = names[name]
    weights = 0
    if function_name == "pt_R2":
        function = pt_R2
    elif function_name == "pt_R":
        function = pt_R
    for num in tqdm(nums, desc=f"Processing Files for {name}", leave=True):
        file = directory + name + "_" + num + ".root"
        
        if type(weights) == type(0):
            weights = function(file, bins)
        else:
            weights += function(file, bins)
        # break

    p0 = [float(weights[0])] + [float(1)]*N #initial guess of parameters
    fitted_beta, covariance = curve_fit(g, x_data, weights, p0=p0,)
    
    # custom_g = lambda x, *args: g(x,*([weights[0]] + list(args)))
    # p0 = [float(1)]*N
    # fitted_beta, covariance = curve_fit(custom_g, x_data,weights,p0=p0)

    plt.scatter(x_data, weights, label='Data')
    y_curve = [g(x, *fitted_beta) for x in x_data]
    plt.plot(x_data,y_curve,'r-', label='Fitted function')
    
    textstr=""
    for j in range(1,5):
        textstr+=r'$M_{%i}=$ %.2e' % (j,moment(j,x_data,weights)) + '\n'
    for j in range(1,len(fitted_beta)):
        textstr+=r'$\beta_{%i}=$ %.2e' % (j,fitted_beta[j]) + '\n'
    
    plt.text(0.7, 0.95, textstr, 
     transform=plt.gca().transAxes, 
     verticalalignment='top', 
     horizontalalignment='left', 
     fontsize=12, 
     bbox=dict(facecolor='white', alpha=0.5))
    plt.ylabel('$P - probability$')
    plt.xlim(0,1)
    
    if function_name == "pt_R":
        plt.xlabel('$R$')
        plt.title("$P_T(R)$ za {}".format(name))
        plt.savefig("./pdfs/pt_R-{}.pdf".format(name))

    elif function_name == "pt_R2":
        plt.xlabel('$R^2$')
        plt.title("$P_T(R^2)$ za {}".format(name))
        plt.savefig("./pdfs/pt_R2-{}.pdf".format(name))
    # plt.show()
    plt.clf()

def g(x, *beta):
    ans = 1.0
    zz = 1.0
    exponent = 0
    try:
        for i in range(1,len(beta)):
            zz *= x
            exponent += (-beta[i]*zz)
    except:
        pass
    ans = np.exp(exponent)
    return ans*beta[0]

#gorozdova 5