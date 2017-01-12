'#financialmath.py'
#%%
import numpy as np

def present_value(cash_flows, time_ids, interest_rates,
                  probs=np.array([1]), power=1):

    for x in [cash_flows, time_ids, interest_rates, probs]:
        if not isinstance(x, np.ndarray):
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
    def __init__(self, x, lx, m=1):
        for x in [cash_flows, time_ids, interest_rates, probs]:
            if not isinstance(x, np.ndarray):
                raise ValueError('Inputs must be numpy arrays.')
        '# extends the age list to include 1/m time steps.'
        k = np.array(np.arange(0,1,1/m)*len(x))
        x = x.repeat(m) + k
        
        '# extends the lives list to include 1/m time steps.'
        dx = lx - np.roll(lx, -1)
        dx[-1] = lx[-1]
        dx = dx.repeat(m)/m
        kk = np.array(np.arange(m)*len(lx))
        lx = lx.repeat(m) - dx*kk
        
        self._x = x
        self._lx = lx
        self._m = m
    
    def omega(self):
        
        return int(self._x[-1] + 1/self._m)
    
    
    def dx(self, x=None):
        d = self._lx - np.roll(self._lx, -1)
        d[-1] = self._lx[-1]
        if x is None:
            return d
        elif (x < 0) | (x > self.omega()*self._m):
            raise ValueError('x must be between 0 and omega.')
        else:
            return d[x]
    
    
    def ptx(self, t, x=None):
        
        if (t < 0) | (t > self.omega()*self._m):
            raise ValueError('t must be between 1 and omega.')
        p = (np.roll(self._lx, -t)/self._lx)
        p[-t:] = 0
        
        if x is None:
            return p
        elif (x < 0) | (x > self.omega()*self._m):
            raise ValueError('x must be between 0 and omega.')
        else:
            return p[x]
    
    def Px(self):
        
        return np.array([self.ptx(i+1) for i in range(self.omega()*self._m)]).T
    
    def qtx(self, t, x=None):
        
        return 1 - self.ptx(t, x)
    
    def Qx(self):
        
        return np.array([self.qtx(i+1) for i in range(self.omega()*self._m)]).T
    
    
    
    def ex(self):
        
        return np.array([[self.ptx(i+1, x) for i in range(self.omega()*self._m)]
                         for x in range(self.omega()*self._m)]).sum(axis=1)
    
    def Lx(self):
        pass
    
    def mtx(self, t, x=None):
        pass
    
    def DataFrame(self):
        
        od = collections.OrderedDict()
        od['x'] = (np.array(self._x))
        od['lx'] = (np.array(self._lx))
        od['px'] = self.ptx(1)
        od['qx'] = self.qtx(1)
        od['dx'] = self.dx()
        od['ex'] = self.ex()
        
        return pd.DataFrame(od, columns=list(od.keys()))
#%%

def Annuity():
    pass