import itertools

def chunks(lst, n, init_pos=0):
    """Yield successive n-sized chunks from lst."""
    for i in range(int(init_pos), len(lst), n):
        yield lst[i:i + n]

def check_spin(ks_file):
   n=0
   result=[]
   if 'KOHN-SHAM MATRIX FOR ALPHA SPIN' in open(ks_file).read():
       n=n+1
   if 'KOHN-SHAM MATRIX FOR BETA SPIN' in open(ks_file).read():
       n=n+1

   if n==2:
      result.append("spin_pol")

   elif n!=0:
      result.append('no_spin_pol')
       
   return result
 
       
def spin_ks_ham(ks_file):
    """
    Function to separate the Hamiltonian based on spin alpha and beta.
    
    Parameter:
    ks_file :: str
      name of the file

    Output:
    result:: tuple
      tuple of lists with the Hamiltonians split up in alpha and beta
      components.
    
    """
    with open(ks_file) as fp, open(ks_file) as fp_1:
        result_alpha=list(itertools.takewhile(lambda x: 'KOHN-SHAM MATRIX FOR BETA SPIN' not in x, 
             itertools.dropwhile(lambda x: 'KOHN-SHAM MATRIX FOR ALPHA SPIN' not in x, fp)))

        result_beta=list(itertools.dropwhile(lambda x: 'KOHN-SHAM MATRIX FOR BETA SPIN' not in x, fp_1))
        

    ks_spin_alpha=list(filter(str.strip,[element.strip('\n') for element in result_alpha]))[1:]

    ks_spin_beta=list(filter(str.strip,[element.strip('\n') for element in result_beta]))[1:]
    
    return ks_spin_alpha, ks_spin_beta

def ks_ham(ks_file):
    """ no spin polarized calculation """
    print ("Hello")
    
# CHECK BETA

#print(check_spin("CH.dat-1_0_78.Log"))

a,b=spin_ks_ham("CH.dat-1_0_78.Log")

# WORKING HAMILTONIAN CHUNKING IN TERMS OF the number of maximum shells +1 (113+1)
ham=list(chunks(a,114))

for i in range(len(ham[0])):
    print (ham[0][i])

