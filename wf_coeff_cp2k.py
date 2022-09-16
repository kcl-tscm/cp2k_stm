import scipy.io
import numpy as np

class wfn_gs(object):
    """ Wave function object, as stored in a .wfn file """

    def __init__(self,
                 natom=None,
                 nspin=None,
                 nao_tot=None,
                 nset_max=None,
                 nshell_max=None):
        #defined at the initialization
        self.natom = int(natom)  # number of atoms
        self.nspin = int(nspin)  # number of spin (i.e., 1 or 2)
        self.nao_tot = int(nao_tot)  # total number of atomic orbitals
        self.nset_max = int(nset_max)  # max number of sets of bs
        self.nshell_max = int(nshell_max)  # mak number of shells in each set

        #spin independent properties
        self.nset = [0 for _i in range(self.natom)]
        self.nshell = [[0 for _i in range(self.nset_max)]
                       for _j in range(self.natom)]
        self.nao = [[[0 for _l in range(self.nshell_max)]
                     for _i in range(self.nset_max)]
                    for _j in range(self.natom)]

        #spin dependent properties
        self.nmo = [0, 0]  # number of molecular orbitals
        self.nocc = [0, 0]  # number of occupied molecular orbitals
        self.nvirt = [0, 0]  # number of virtual molecular orbitals
        self.nel = [0, 0]  # number of electrons

    def initialize_lists(self, ispin):
        # Initialize (or reset to zero) alpha or beta molecular orbitals
        if ispin == 0:
            # eigenvalues for each molecular orbital
            self.eigen = [[0 for _i in range(self.nmo[ispin])]]
            # occupancies for each molecular orbital
            self.occup = [[0 for _i in range(self.nmo[ispin])]]
            #coefficient list for each moleculat orbital
            self.coeff = [[[0 for _j in range(self.nao_tot)]
                           for _i in range(self.nmo[ispin])]]
        elif ispin == 1:
            self.eigen.append([0 for _i in range(self.nmo[ispin])])
            self.occup.append([0 for _i in range(self.nmo[ispin])])
            self.coeff.append([[0 for _j in range(self.nao_tot)]
                               for _i in range(self.nmo[ispin])])

    def makeopenshell(self):
        # Makes a restricted (e.g., nspin=1) wfn open shell (e.g., nspin=2)
        if self.nspin == 2:
            return
        elif self.nspin == 1:
            self.nspin = 2
            self.nel[0] /= 2
            self.eigen[0] = [x / 2. for x in self.eigen[0]]
            self.occup[0] = [x / 2. for x in self.occup[0]]
            self.nmo[1] = self.nmo[0]
            self.nocc[1] = self.nocc[0]
            self.nvirt[1] = self.nvirt[0]
            self.nel[1] = self.nel[0]
            self.initialize_lists(1)
            self.eigen[1] = self.eigen[0]
            self.occup[1] = self.occup[0]
            self.coeff[1] = self.coeff[0]
            return

    def add_nel(self, nel_geom, charge, multiplicity):
        #default: charge=0, mult=1
        # Check that nel/mult/charge are consistent
        # If nspin=2 split into alpha and beta
        # If charge!=0, add/remove electrons
        # If multiplicity>1, rearrange alpha and beta electrons
        self.charge = charge
        self.mult = multiplicity
        if self.nspin == 1:
            if self.mult > 1:
                print("WARNING: nspin=1 but multiplicity>1! EXIT")
                sys.exit()
            self.nel[0] = int(nel_geom - self.charge)
            if self.nel[0] % 2 == 1:
                print("WARNING: nspin=1 but odd number of electrons! EXIT")
                sys.exit()
            self.nmo[0] = int(self.nel[0] / 2.0)  # 2 electr. per MO
            self.nocc[0] = self.nmo[0]  # all the orbitals are occupied
            self.nvirt[0] = 0  # no virtual orbitals
            self.initialize_lists(0)
            for i in range(self.nmo[0]):
                self.eigen[0][i] = 0.0
                self.occup[0][i] = 2.0
        elif self.nspin == 2:
            if (nel_geom - self.charge % 2 == 1) and (self.mult % 2 == 1):
                print("WARNING: both number of electrons and multiplicity \
                       are odd! EXIT")
            self.nel[0] = int((nel_geom - self.charge) / 2.0 +
                              (self.mult - 1) / 2.0)
            self.nel[1] = int((nel_geom - self.charge) / 2.0 -
                              (self.mult - 1) / 2.0)
            for ispin in range(self.nspin):
                self.nmo[ispin] = self.nel[ispin]  # 1 electr. per MO
                self.nocc[ispin] = self.nmo[
                    ispin]  # all the orbitals are occup.
                self.nvirt[ispin] = 0  # no virtual orbitals
                self.initialize_lists(ispin)
                for i in range(self.nmo[ispin]):
                    self.eigen[ispin][i] = 0.0
                    self.occup[ispin][i] = 1.0

def read_wfn_gs_file(wfn_file):
    """ Parser for the .wfn file, that returns the wfn object """
    inpfile = scipy.io.FortranFile(wfn_file, "r")
    natom, nspin, nao_tot, nset_max, nshell_max = inpfile.read_ints()
    w = wfn_gs(natom, nspin, nao_tot, nset_max, nshell_max)
    nset = inpfile.read_ints()
    for i, val in enumerate(nset):
        w.nset[i] = val
    nshell = inpfile.read_ints()
    for i, val in enumerate(nshell):
        iatom, iset = divmod(i, w.nset_max)
        w.nshell[iatom][iset] = val
    nao = inpfile.read_ints()
    for i, val in enumerate(nao):
        iatom, isetshell = divmod(
            i, w.nset_max * w.nshell_max
        )  # this splits the N of atomic orbitals in one list for each atom
        iset, ishell = divmod(
            isetshell, w.nshell_max
        )  # this splits the N of atomic orbitals in one list for each set
        w.nao[iatom][iset][ishell] = val
    for ispin in range(w.nspin):
        nmo, nocc, nvirt, nel = inpfile.read_ints()
        w.nmo[ispin] = nmo
        w.nocc[ispin] = nocc
        w.nvirt[ispin] = nvirt
        w.nel[ispin] = nel
        w.initialize_lists(ispin)
        eigen_and_occup = inpfile.read_reals()  #this line contains both
        for i in range(w.nmo[ispin]):
            w.eigen[ispin][i] = eigen_and_occup[i]
        for i in range(w.nmo[ispin]):
            w.occup[ispin][i] = eigen_and_occup[i + w.nmo[ispin]]
        for imo in range(w.nmo[ispin]):
            coeff = inpfile.read_reals()
            for i, val in enumerate(coeff):
                w.coeff[ispin][imo][i] = val
    inpfile.close()
    return w

#def read_wfn_exc_file(wfn_file):
#    inpfile = scipy.io.FortranFile(wfn_file, "r")
#    num_exc_sta, nspin, nao_tot, nshell_max = inpfile.read_ints()
#    w = wfn(natom, nspin, nao_tot, nset_max, nshell_max)
    


if __name__ == "__main__":
    total=read_wfn_gs_file("BENZENE-UNPERTURBED-RESTART.wfn.bak-3")

    # Getting number of number of molecular orbitals
    #up, down=total.nmo

    print (total.nao_tot)
    #Printing coefficients
    #for idx, values in enumerate(total.coeff[0]):
    #    print (idx,len(values))
    
    # Printing occupations
    #print (total.occup)

    # Printing eigenvalues
    #print (total.eigen)
