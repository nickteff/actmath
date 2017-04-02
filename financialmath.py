'#financialmath.py'
#%%
import numpy as np

def present_value(cash_flows, time_ids, interest_rates,
                  probs=np.array([1]), power=1):

    for item in [cash_flows, time_ids, interest_rates, probs]:
        if not isinstance(item, np.ndarray):
            raise ValueError('Inputs must be numpy arrays.')

    if len(cash_flows) != len(time_ids):
        raise ValueError('Cash flow and time ids are not the same length.')

    if len(interest_rates) == 1:
        ir = interest_rates[0]
        interest_rates = np.ones_like(time_ids)*ir

    elif len(probs) == 1:
        p = probs[0]
        probs = np.ones_like(cash_flows)*p

    elif len(cash_flows) != len(probs):
        raise ValueError('Cash flow and probabilities are '
                         'not the same length.')

    '#discount values'
    v = (1 + interest_rates)**(time_ids)
    v = 1/v

    return (np.power(cash_flows, power)*np.power(v, power)*probs).sum()

#%%
import collections
import numpy as np
import pandas as pd


class Lifetable(object):
    def __init__(self, x, lx, k=1):
        for item in [x, lx]:
            if not isinstance(item, np.ndarray):
                raise ValueError('Inputs must be numpy arrays.')
                
        '# extends the age list to include 1/k time steps.'
        kk = np.array(list(np.arange(0, 1, 1/k))*len(x))
        x = x.repeat(k) + kk
        
        '# extends the lives list to include 1/k time steps.'
        dx = lx - np.roll(lx, -1)
        dx[-1] = lx[-1]
        dx = dx.repeat(k)/k
        kk = np.array(list(np.arange(k))*len(lx))
        lx = lx.repeat(k) - dx*kk
        
        self._x = x
        self._lx = lx
        self._k = k
    
    def omega(self):
        
        return int(self._x[-1] + 1/self._k)

    def dx(self, x=None):
        d = self._lx - np.roll(self._lx, -1)
        d[-1] = self._lx[-1]
        if x is None:
            return d
        elif (x < 0) | (x > self.omega()*self._k):
            raise ValueError('x must be between 0 and omega.')
        else:
            return d[x]

    def ptx(self, t, x=None):

        if (t < 0) | (t > self.omega()*self._k):
            raise ValueError('t must be between 1 and omega.')
        p = (np.roll(self._lx, -t)/self._lx)
        p[-t:] = 0

        if x is None:
            return p
        elif (x < 0) | (x > self.omega()*self._k):
            raise ValueError('x must be between 0 and omega.')
        else:
            return p[int(x*self._k)]

    def Px(self):
        
        return np.array([self.ptx(i+1) for i in range(self.omega()*self._k)]).T

    def qtx(self, t, x=None):
        
        return 1 - self.ptx(t, x)

    def Qx(self):
        
        return np.array([self.qtx(i+1) for i in range(self.omega()*self._k)]).T

    def ex(self):
        
        return np.array([self.ptx(i+1)
                         for i in np.arange(self.omega()*self._k)]).sum(axis=0)
    
    def DataFrame(self):
        
        od = collections.OrderedDict()
        od['x'] = (np.array(self._x))
        od['lx'] = (np.array(self._lx))
        od['px'] = self.ptx(1)
        od['qx'] = self.qtx(1)
        od['dx'] = self.dx()
        od['ex'] = self.ex()

        return pd.DataFrame(od, columns=list(od.keys()), index=od['x'])
#%%

def annuity(n, interest_rates, m, age=None, lifetable=None, 
            k=1, power=1, payment='due'):
    #include error checking on age, n, m, k, payment 
    #make the default payment immediate
    #computation of quantities, assuming fractional payments
    if (age is not None) & (lifetable is not None):
        if k != lifetable._k:
            raise ValueError('The fractional timestep k is not the same as the '
                             'lifetable.')
            
        n = np.min([(lifetable.omega()-age-m)*k, n*k])
        probs = lifetable.ptx(1)[age:age+n]
    else:
        n = n*k
        probs = np.array([1])
    
    cash_flows = np.array(1/k).repeat(n)
    time_ids = np.arange(0, n/k, 1/k) + m 

    if payment == "immediate":
        time_ids = time_ids + 1/k

    return present_value(cash_flows=cash_flows,
                         time_ids=time_ids,
                         interest_rates=interest_rates,
                         probs=probs,
                         power=power)
