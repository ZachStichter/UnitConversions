from typing import Mapping

FORCE_MODE = None

def force_units(forcetype:str|None=None)->str|None:
    '''
    Set a force mode; that is, for numbers that automatically inherit their unit system, force the unit system to be either SI or AU

    Args:
        forcetype (str|None): the type of system to force. This can be either "SI" (SI units), "AU" (atomic units), or None (inherit dynamically)
    '''
    global FORCE_MODE
    if forcetype == None:
        FORCE_MODE = None
        return FORCE_MODE
    if forcetype.upper().strip() in ['SI', 'AU']:
        FORCE_MODE = forcetype.upper().strip()
    else:
        FORCE_MODE = None
    return FORCE_MODE

class DimensionalQuantity(float):
    '''
    An Extension of float to allow direct conversion with dimensional units.

    Args:
        value (float): the magnitude of the dimensional quantity to store
        dims (dict): a dictionary of the form {'unit':quantity}, where the units associated with value are unit^quantity.
            Example: 3.14 m/s -> DimensionalQuantity(3.14, {'m':1,'s':-1})
    '''
    def __new__(cls, value: float, dims: Mapping[str, int]):
        obj = super().__new__(cls, value)
        obj._dims = dict(dims)
        return obj

    def to_AU(self) -> float:
        '''
        Force convert the dimensional quantity to atomic units (hbar=1, etc.) from any valid SI or extended units
        '''
        return SI.to_AU(self, self._dims)

    def from_AU(self) -> float:
        '''
        Force convert the dimensional quantity from atomic units (hbar=1, etc.) to SI base units
        '''
        return SI.from_AU(self, self._dims)

    def to_Force(self)->float:
        '''
        Using the FORCE_MODE parameter (set at startup),
        '''
        global FORCE_MODE
        match FORCE_MODE:
            case 'AU':
                return SI.to_AU(self, self._dims)
            case _:
                return SI.to_Base(self, self._dims)

    def to_Base(self)->float:
        return SI.to_Base(self, self._dims)

    def from_Molar(self)->float:
        return SI.from_Molar(self, self._dims)

    @property
    def AU(self):
        return self.to_AU()

    @property
    def base(self):
        return self.to_Base()

    @property
    def forced(self):
        return self.to_Force()

    @property
    def dims(self)->Mapping:
        return self._dims

    @property
    def fdims(self) -> str:
        fstr = ''
        for k, v in {k: v for k, v in sorted(self._dims.items(), key=lambda item: item[1], reverse=True)}.items():
            fstr += f'*{k}' if v > 0 else f'/{k}'
            if abs(v) > 1:
                fstr += f'^{abs(v)}'
        return fstr[1:]

    def float(self):
        return float(self)

    @property
    def flt(self) -> float:
        return self.float()

    def __mul__(self, other):
        if isinstance(other, DimensionalQuantity):
            new_val = float(self) * float(other)
            new_dims = self._combine_dims(self._dims, other._dims, op='+')
            return DimensionalQuantity(new_val, new_dims)
        else:
            return DimensionalQuantity(float(self) * other, self._dims)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, DimensionalQuantity):
            new_val = float(self) / float(other)
            new_dims = self._combine_dims(self._dims, other._dims, op='-')
            return DimensionalQuantity(new_val, new_dims)
        else:
            return DimensionalQuantity(float(self) / other, self._dims)

    def __rtruediv__(self, other):
        if isinstance(other, DimensionalQuantity):
            raise NotImplementedError("DimensionalQuantity / DimensionalQuantity already handled")
        new_val = other / float(self)
        new_dims = {k: -v for k, v in self._dims.items()}
        return DimensionalQuantity(new_val, new_dims)

    def __pow__(self, power):
        if isinstance(power, DimensionalQuantity):
            if power._dims:  # Non-dimensionless exponent
                raise ValueError("Exponent must be dimensionless")
            power = float(power)

        if not isinstance(power, (int, float)):
            return NotImplemented

        new_val = float(self) ** power
        new_dims = {k: v * power for k, v in self._dims.items()}
        return DimensionalQuantity(new_val, new_dims)

    def __rpow__(self, base):
        # base ** self
        if self._dims:
            raise ValueError("Exponent must be dimensionless")
        return base ** float(self)


    def _combine_dims(self, d1, d2, op='+'):
        result = d1.copy()
        for k, v in d2.items():
            if k in result:
                result[k] += v if op == '+' else -v
                if result[k] == 0:
                    del result[k]
            else:
                result[k] = v if op == '+' else -v
        return result

class _SI():
    '''
    USING UPDATED CODATA MATERIAL 2022

    REF: 10.1103/RevModPhys.97.025002
    '''
    pi =            3.1415926535897932384626433832
    rad =           DimensionalQuantity(1, {'rad':1})
    sr =            DimensionalQuantity(1, {'sr':1})
    h =             DimensionalQuantity(6.62607015E-34, {'kg':1, 'm':2, 's':-1})
    e =             DimensionalQuantity(1.602176634E-19, {'A':1, 's':1})
    k_b =           DimensionalQuantity(1.380649E-23, {'kg':1, 'm':2, 'K':-1, 's':-2})
    N_A =           DimensionalQuantity(6.02214076E23, {'mol':-1})
    c =             DimensionalQuantity(299792458, {'m':1, 's':-1})
    v_Cs =          DimensionalQuantity(9192631770, {'s':-1})
    K_cd =          DimensionalQuantity(683, {'cd':1, 'sr':1, 's':3, 'kg':-1, 'm':-2})
    s =             float(v_Cs)/v_Cs
    m =             (float(v_Cs)*c)/(float(c)*v_Cs)
    kg =            (float(c*c))*(h*v_Cs)/(float(h)*float(v_Cs)*c*c)
    A =             (e*v_Cs)/(float(e)*float(v_Cs))
    C =             (e)/(float(e))
    K =             (float(k_b))/float(h*v_Cs)*(h*v_Cs)/k_b
    mol =           float(N_A)/N_A
    cd =            (K_cd*h*v_Cs*v_Cs/sr)/float(K_cd*h*v_Cs*v_Cs)
    
class _AMU(_SI):
    hbar = _SI.h / (2*_SI.pi)
    m_e = _SI.kg*9.1093837139E-31
    mu_0 = _SI.pi*_SI.kg*_SI.m/_SI.s/_SI.s/_SI.A/_SI.A*0.99999999987E-7*4
    eps_0 = 1/_SI.c/_SI.c/mu_0
    k_0 = _SI.pi*eps_0*4
    a_0 = k_0*hbar*hbar/m_e/_SI.e/_SI.e
    alpha = _SI.e*_SI.e/k_0/hbar/_SI.c
    E_h = alpha*alpha*m_e*_SI.c*_SI.c
    t_0 = hbar/E_h
    _conv_AMU_T = DimensionalQuantity(1/t_0, {'AMU.T':1, 's':-1})
    _conv_AMU_L = DimensionalQuantity(1/a_0, {'AMU.L':1, 'm':-1})
    _conv_AMU_M = DimensionalQuantity(1/m_e, {'AMU.M':1, 'kg':-1})
    _conv_AMU_W = DimensionalQuantity(1/E_h, {'AMU.W':1, 'J':-1})
    _conv_AMU_Q = DimensionalQuantity(1/_SI.e, {'AMU.Q':1, 'C':-1})
    _conv_AMU_N = DimensionalQuantity(_SI.N_A, {'AMU.N':1, 'mol':-1})
    _conv_eV_J = DimensionalQuantity(float(_SI.e), {'J':1, 'eV':-1})


    @classmethod
    def to_AU(cls, quant:float, qdims:Mapping[str, int])->DimensionalQuantity:
        '''
        Given a quantity in any SI units, convert it to base SI, then to atomic units.

        Select from the following list:
            mass (g)
            length (m)
            time (s)
            charge (C)
            energy (J|eV)
            number (mol)

        Args:
            quant: the quantity to convert, in valid SI units
            quant_units: a dictionary containing mappings of the form {unit:exponent},
                where unit is the fundamental class, and exponent is the current exponent. e.g.: m/s->{"m":1, "s":-1}

        Returns:
            A DimensionalQuantity with magnitude of quant (converted to atomic units) and qdims of the corresponding atomic unit.
        '''
        quant = DimensionalQuantity(quant, qdims).base
        modified_unit = quant
        dims = quant.dims
        for key, val in dims.items():
            conv_factor = 1
            match key.lower():
                case 'kg':
                    conv_factor = _AMU._conv_AMU_M
                case 'm':
                    conv_factor = _AMU._conv_AMU_L
                case 's':
                    conv_factor = _AMU._conv_AMU_T
                case 'c':
                    conv_factor = _AMU._conv_AMU_Q
                case 'ev':
                    conv_factor = _AMU._conv_eV_J
                    conv_factor *= _AMU._conv_AMU_W
                case 'j':
                    conv_factor = _AMU._conv_AMU_W
                case 'mol':
                    conv_factor = _AMU._conv_AMU_N
                case 'a':
                    conv_factor = DimensionalQuantity(_AMU._conv_AMU_Q/_AMU._conv_AMU_T,{'A':-1, 'AMU.Q':1, 'AMU.T':-1})
            modified_unit = modified_unit * conv_factor**(val)
        return modified_unit

    @classmethod
    def from_AU(cls, quant:float, qdims:Mapping[str, int])->DimensionalQuantity:
        '''
        Given a quantity in base atomic units, convert it to base SI.

        Select from the following list:
            mass (AMU.M)
            length (AMU.L)
            time (AMU.T)
            charge (AMU.Q)
            energy (AMU.W)
            number (AMU.N)

        Args:
            quant: the quantity to convert, in valid SI units
            quant_units: a dictionary containing mappings of the form {unit:exponent},
                where unit is the fundamental class, and exponent is the current exponent. e.g.: m/s->{"m":1, "s":-1}

        Returns:
            A DimensionalQuantity with magnitude of quant (converted to atomic units) and qdims of the corresponding atomic unit.
        '''
        mod = DimensionalQuantity(quant, qdims)
        for key, val in qdims.items():
            conv_factor = 1
            match key.lower():
                case 'amu.m':
                    conv_factor = 1/_AMU._conv_AMU_M
                case 'amu.l':
                    conv_factor = 1/_AMU._conv_AMU_L
                case 'amu.t':
                    conv_factor = 1/_AMU._conv_AMU_T
                case 'amu.q':
                    conv_factor = 1/_AMU._conv_AMU_Q
                case 'amu.w':
                    conv_factor = 1/_AMU._conv_AMU_W
                case 'amu.n':
                    conv_factor = 1/_AMU._conv_AMU_N
            mod *= conv_factor**(val)
        return mod

class SI(_AMU):
    @classmethod
    def parse(cls, unit:str, sys:str='f')->DimensionalQuantity:
        '''
        Given a string of form "num unit" return the float representing that unit in base SI or A.U.

        Args:
            unit (str): A string of the form "num unit" where num is the scalar quantity and unit is the SI unit
            sys (str): The desired unit system for the output: f|a|b
                f: force; use the forced units
                a: AU; use atomic units
                b: base; use SI base units
        '''
        USE_FORCE = False
        USE_BASE = True
        USE_AU = False
        match sys.lower():
            case 'f':
                USE_FORCE = True
                USE_BASE = False
                USE_AU = False
            case 'a':
                USE_FORCE = False
                USE_BASE = False
                USE_AU = True
            case _:
                USE_FORCE = False
                USE_AU = False
                USE_BASE = True
        val, dim = unit.split(maxsplit=1)
        if ' ' in dim.strip():
            dim = dim.split()
        else:
            dim = [dim]
        for idx, exp in enumerate(dim):
            dim[idx] = exp.split('^')
            if len(dim[idx]) == 1:
                dim[idx].append(1)
            dim[idx][1] = int(dim[idx][1])
        if USE_AU:
            return DimensionalQuantity(val, {un:exp for un, exp in dim}).AU
        elif USE_FORCE:
            return DimensionalQuantity(val, {un:exp for un, exp in dim}).forced
        else:
            return DimensionalQuantity(val, {un:exp for un, exp in dim}).base

    @classmethod
    def to_Base(cls, val, dims):
        q = DimensionalQuantity(val, dims).from_AU()
        val = q.flt
        dims = q.dims
        fundamental_units = ['s', 'm', 'g', 'A', 'K', 'mol', 'cd']
        conversions = {'eV':[_SI.e.flt,'J']}
        extended_units = list(conversions.keys())
        prefixes = {'a':-18, 'f':-15, 'p':-12, 'n':-9, 'u':-6, 'm':-3,'c':-2,'d':-1,'da':1,'h':2,'k':3, 'M':6,'G':9,'T':12,'P':15,'E':18}
        new_dims = {}
        for k, v in dims.items():
            conv_factor = 0
            prefactor = 1
            q = ''
            u_ = ''
            for u in fundamental_units+extended_units:
                if k.endswith(u):
                    q = k.split(u)[0]
                    u_ = u
                    break
            if u_ in extended_units:
                prefactor *= conversions[u_][0]
            if q != '':
                try:
                    assert q in prefixes
                except AssertionError:
                    raise ValueError(f'{q} is not a valid SI prefix. If you have used the greek mu, use u instead.')
                if u_ in extended_units:
                    prefactor *= 10**prefixes[q]
                else:
                    conv_factor += prefixes[q]*v
            else:
                new_dims[k] = v
            if u_ == 'g':
                conv_factor -= 3*v
                new_dims['kg'] = v
            elif u_ in extended_units:
                new_dims[conversions[u_][1]] = v
                new_dims.pop(k)
            elif u_ != '':
                new_dims[u_] = v
        return DimensionalQuantity(prefactor*val*10**conv_factor,new_dims)

    @classmethod
    def from_Molar(cls, val, dims):
        if 'mol' in dims.keys():
            return DimensionalQuantity(1,{'mol':dims['mol']}).AU*val
        else:
            return DimensionalQuantity(val, dims)

if __name__ == '__main__':
    print(SI.s, SI.s.fdims)
    print(SI.m, SI.m.fdims)
    print(SI.kg, SI.kg.fdims)
    print(SI.A, SI.A.fdims)
    print(SI.C, SI.C.fdims)
    print(SI.K, SI.K.fdims)
    print(SI.mol, SI.mol.fdims)
    print(SI.cd, SI.cd.fdims)
    print(SI.hbar, SI.hbar.fdims)
    print(SI.mu_0, SI.mu_0.fdims)
    print(SI.eps_0, SI.eps_0.fdims)
    print(SI.k_0, SI.k_0.fdims)
    print(SI.a_0, SI.a_0.fdims)
    print(SI.alpha, SI.alpha.fdims)
    print(SI.E_h, SI.E_h.fdims)
    print(SI.t_0, SI.t_0.fdims)
    print(SI._conv_AMU_T, SI._conv_AMU_T.fdims)
    a = DimensionalQuantity(1,{'nm':1})
    print(a.AU, a.AU.fdims)
    b = DimensionalQuantity(1,{'kg':1,'m':2,'s':-2})
    print(1/b.AU, b.AU.fdims)
    print(_AMU.hbar.AU, _AMU.hbar.AU.fdims)
    print(_AMU.m_e.AU, _AMU.m_e.AU.fdims)
    print(_AMU.t_0.AU, _AMU.t_0.AU.fdims)
    print(_AMU.k_0.AU, _AMU.k_0.AU.fdims)
    c = DimensionalQuantity(1,{'eV':1})
    print(c.AU, c.AU.fdims)
    d = SI.parse('16.46 eV', 'b')
    print(d)
    print(d.flt)