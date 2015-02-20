# Tools for working with Stanley sequences
# - Calculating terms
# - Reading and writing in base 3 (see paper A)
# - Independent sequences (see papers A and B)
# - Modular sequences (see paper C)
# - Basic sequences (see paper C)
#
# Papers:
#
# A) D. Rolnick, "On the classification of Stanley sequences," preprint arXiv:1408.1940.
# B) D. Rolnick and P. S. Venkataramana, "On the growth of Stanley sequences," preprint arXiv:1408.4710.
# C) R. Moy and D. Rolnick, "Novel structures in Stanley sequences," in preparation.

from numpy import base_repr,log
from fractions import gcd

class StanSeq:
    #gen_set should be a list of increasing integers, starting with 0
    def __init__(self,gen_set,base_3=False):
        assert type(gen_set)==list and len(gen_set)!=0 and gen_set[0]==0, "Generating set should be a list of increasing integers, starting with 0."
        assert test_3_free(gen_set), "Generating set is not 3-free."
            
        if base_3:
            A=from_base_3(gen_set)
        else:
            A=gen_set
    
        #terms and covered are stored as bits within a Python long integer
        #a 1 in covered means that a number is *ruled out*
        #terms is stored with leftmost bit representing 0, rightmost n
        #covered is stored with leftmost bit representing 2n, rightmost n+1
        #thanks to Ricky Liu for this clever idea
        
        terms=1
        covered=0
        next_in_set=1
        for current in range(1,A[-1]+1):
            if current==A[next_in_set]:
                terms=(terms<<1)+1
                covered=(terms|covered)>>1    
                next_in_set+=1
            else:
                terms<<=1
                covered>>=1
        
        self.terms=terms
        self.covered=covered
        self.max_tested=A[-1]
        self.num_terms=next_in_set-1
        self.gen_set=A
    
    #returns the sequence as a list
    def as_list(self):
        terms=self.terms
        current=self.max_tested
        the_list=[]
        while terms!=0:
            if terms & 1:
                the_list.insert(0,current)
            terms>>=1
            current-=1
        return the_list
    
    #tests the next integer for admission to the sequence and adds if possible
    def extend(self):
        if self.covered & 1:
            self.terms<<=1
            self.covered>>=1
            self.max_tested+=1
        else:
            self.terms=(self.terms<<1)+1
            self.covered=(self.terms|self.covered)>>1
            self.num_terms+=1
            self.max_tested+=1
        return self    
                                        
    #returns the first num_terms terms of the sequence
    #if options = "list (base 10)", then output is a list in base 10
    #if options = "list (base 3)", then output is a list in base 3
    #if options = "StanSeq", then output is a StanSeq object
    def n_terms(self,num_terms,options="list (base 10)"):
        while self.num_terms<num_terms:
            self.extend()
        if options=="list (base 10)":
            return self.as_list()[:num_terms+1]
        if options=="list (base 3)":
            return to_base_3(self.as_list()[:num_terms+1])
        if options=="StanSeq":
            return self
    
    #returns all terms of the sequence up to max_term
    #if options = "list (base 10)", then output is a list in base 10
    #if options = "list (base 3)", then output is a list in base 3
    #if options = "StanSeq", then output is a StanSeq object
    def up_to_n(self,max_term,options="list (base 10)"):
        while self.max_tested<max_term:
            self.extend()
        if options=="list (base 10)":
            return [x for x in self.as_list() if x<=max_term]
        if options=="list (base 3)":
            return to_base_3([x for x in self.as_list() if x<=max_term])
        if options=="StanSeq":
            return self
        
    #returns the minimal generating set A for the sequence
    def min_gen_set(self):
        #tests if cutting off the last element preserves the sequence
        while self.gen_set!=[0] and StanSeq(self.gen_set[:-1]).n_terms(len(self.gen_set)-1)==self.gen_set:
            self.gen_set=self.gen_set[:-1]
        return self.gen_set
        
    #tests whether the sequence is independent up to n terms, and if so optionally returns parameters
    def test_indep(self,num_terms,parameters=True):
        indep=False
        max_power=int(log(num_terms)/log(2))
        power=0
        while 2**power<len(self.min_gen_set()):
            power+=1
        while power<=max_power and not indep:
            S=self.n_terms(2**power)
            indep=test_mod_set(S[:-1],S[-1])
            power+=1
        if indep and parameters:
            char=2*S[-2]-S[-1]+1
            scal_numer=S[-1]/gcd(S[-1],(3**(power-1)))
            scal_denom=(3**(power-1))/gcd(S[-1],(3**(power-1)))
            return "Independent with character %d and scaling factor %s." %(char,(str(scal_numer)+'/'+str(scal_denom)) if scal_denom!=1 else str(scal_numer))
        else:
            return indep
    
    #tests whether the sequence is modular up to n terms, and if so optionally returns parameters
    def test_modular(self,num_terms,parameters=True):
        modular=False
        min_index=len(self.min_gen_set())
        index=min_index
        while index<=num_terms and not modular:
            S=self.n_terms(index)
            if test_mod_set(S[:-1],S[-1])==True:
                modular=True
            index+=1
        if modular and parameters:
            return "Modular with modulus %d and modular set %s." %(S[-1],'['+','.join(str(x) for x in S[:-1])+']')
        else:
            return modular
        
    #tests whether the sequence is basic up to n terms, and if so optionally returns basis
    def test_basic(self,num_terms,parameters=True):
        basic=True
        trial_basis=[]
        S=self.n_terms(num_terms)
        index=1
        while index<=num_terms and basic:
            predictions=set_of_sums(trial_basis)
            if len(predictions)<=index or predictions[index]>S[index]:
                trial_basis.append(S[index])
            elif predictions[index]<S[index]:
                basic=False
            index+=1
        if basic and parameters:
            return "Basic with basis %s." %('['+','.join(str(x) for x in trial_basis) +',...]')
        else:
            return basic

#inputs a list of base 10 integers and outputs a list of integers with digits 0,1,2
def to_base_3(alist):
    output=[]
    for base_3 in alist:
        output.append(int(base_repr(base_3,3)))
    return output
    
#inputs a list of integers with digits 0,1,2 and outputs a list of base 10 integers
def from_base_3(alist):
    output=[]
    for base_3 in alist:
        base_10=0
        for digit in str(base_3):
            base_10=3*base_10+int(digit)
        output.append(base_10)
    return output

#tests if a list is 3-free
def test_3_free(alist):
    alist=sorted(alist)
    length=len(alist)
    for i in range(length):
        for j in range(i+1,length):
            if 2*alist[j]-alist[i] in alist:
                return False
    return True

#tests if a list is a modular set for a given modulus
#all elements of alist must be smaller than modulus and in increasing order
def test_mod_set(alist,modulus):
    length=len(alist)
    covered=0
    for x in alist:
        covered |= (1<<x)
    for i in range(length):
        for j in range(i+1,length):
            x=(2*alist[j]-alist[i]) % modulus
            covered |= (1<<x)
            if x in alist:
                return "Set is not 3-free modulo %d." %(modulus)
    if covered+1==(1<<modulus):
        return True
    else:
        return False

#constructs the (ordered) multiset of sums of subsets of the given list
def set_of_sums(alist):
    length=len(alist)
    output=[0]
    for x in alist:
        output_copy=output[:]
        for y in output_copy:
            output.append(x+y)
    return sorted(output)

#tests if a list is a valid basis
def test_basis(alist):
    if test_mod_set(set_of_sums(alist),(alist[-1])*3)==True:
        return True
    else:
        return False

#constructs a sequence with a given basis
#later terms of basis inferred: e.g. [7,6,9] read as (7,6,9,27,81,243,729,...) 
def basic_seq(basis):
    if test_basis(basis):
        return StanSeq(set_of_sums(basis))
    else:
        return "Invalid basis"
        
#constructs a modular set that this is the product of two modular sets
#the new set is modular with respect to the product of the original moduli
def prod_mod_sets(alist1,modulus1,alist2,modulus2):
    assert test_mod_set(alist1,modulus1)==True and test_mod_set(alist2,modulus2)==True, "Sets must be modular with respect to given moduli."
    output=[]
    for x in alist1:
        for y in alist2:
            output.append(x+modulus1*y)
    return sorted(output)


## TESTS
#print
#print "Construct the terms up to 100 of the Stanley sequence S(0):"
#print StanSeq([0]).up_to_n(100)
#print
#print "Construct the first 100 terms of the Stanley sequence S(0):"
#print StanSeq([0]).n_terms(100)
#print
#print "Construct the terms up to 81 of the Stanley sequence S(0,1,100), where input is given in ternary:"
#print StanSeq([0,1,100],True).up_to_n(81)
#print
#print "Construct the terms up to 81 of the Stanley sequence S(0,1,9), where output is given in ternary:"
#print StanSeq([0]).up_to_n(81,"list (base 3)")
#print
#print "Return the minimal generating set of the Stanley sequence S(0,1,7,8,10,11,17):"
#print StanSeq([0,1,7,8,10,11,17]).min_gen_set()
#print
#print "Test whether the Stanley sequence S(0,1,7) is independent (using 100 terms), and return parameters:"
#print StanSeq([0,1,7]).test_indep(100)
#print
#print "Test whether the Stanley sequence S(0,1,7) is modular (using 100 terms), and return parameters:"
#print StanSeq([0,1,7]).test_modular(100)
#print
#print "Test whether the Stanley sequence S(0,1,7) is basic (using 100 terms), and return basis:"
#print StanSeq([0,1,7]).test_basic(100)
#print
#print "Test whether (7,6,9,27,...) is a valid basis:"
#print test_basis([7,6,9,27])
#print
#print "Construct 100 terms of a Stanley sequence with basis (7,6,9,27,...):"
#print basic_seq([7,6,9,27]).n_terms(100)