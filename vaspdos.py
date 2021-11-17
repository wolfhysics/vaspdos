# code by Farrel Dzaudan Naufal, 13318048
# CMD-QE Laboratory, Engineering Physics Department, Institut Teknologi Bandung
# for free use and distribution
# Rev 17/11/2021, Ver 1.0       Currently only works for:
#                               1. Plotting 1 atom LDOS for s and p orbitals 
#                               2. Calculation must be spin polarized


def main(argv = None):
    # initiate some variables
    ion_num = 0
    spin_component = 0
    kpts_total = 0
    bands_total = 0
    ions_total = 0
    pdos_car = False
    plot_data = True
    plot_lmax = 1
    smearing = 0.1
    lm_orbitals = None
    ion_pdos_data = []
    # set default of argv to None
    # if this module is imported to another module, input goes in argv variable
    # if module run by terminal, main() is executed, and argv = None, so read the terminal args by sys.argv
    if argv is None:
        argv = sys.argv
    if len(argv) < 2:
        return "No argument given, terminating program."
    PROCAR_obj = open('PROCAR', 'r')
    PROCAR_lines = PROCAR_obj.readlines()
    PROCAR_obj.close()

    # read argvs
    try:
        for i in range(len(argv)-1):
            if argv[i] == '-ion':
                if argv[i+1][0].upper() == 'A':
                    ion_num = 0
                    return "This feature is not supported yet, terminating program."
                else:
                    ion_num = int(argv[i+1])
                    file_out_name = 'ION%d.pdos'%ion_num
            if argv[i] == '-smearing':
                smearing = float(argv[i+1])
            if argv[i] == '-lmax':
                plot_lmax = int(argv[i+1])
        print("Generating PDOS plot of")
        print("    ion no %d"%ion_num)
        print("    lmax = %d"%plot_lmax)
        print("    smearing = %f"%smearing)
    except:
        return "Error in reading arguments, terminating program."

    # Read all lines in PROCAR file
    for line in PROCAR_lines:
        if len(line) > 2:

            # get number of k-points, bands, and ions
            if '# of k-points:' in line:
                spin_component += 1
                kpts_total = int(line.split()[3])
                bands_total = int(line.split()[7])
                ions_total = int(line.split()[-1])
                print('Spin component = %d\ntotal bands = %d\ntotal ions = %d'%(spin_component,bands_total,ions_total))
            
            # get band number, energy, and occupation
            if 'band' and 'energy' and 'occ.' in line:
                band_num = int(line.split()[1])
                band_en = float(line.split()[4])
                band_occ = int(float(line.split()[-1]))
            
            # read lm-decomposed PDOS
            if line.split()[0] == 'ion':
                pdos_car = True
                pdos_matrix = []
                if lm_orbitals == None:
                    lm_orbitals = line.split()[1:-1]
                    lm_orb_tot = len(lm_orbitals)

            if pdos_car:
                pdos_matrix.append(line.split())
            
            if line.split()[0] == 'tot':
                pdos_car = False
                norm_const = float(pdos_matrix[-1][-1])
                for j in range(1,ions_total+1):
                    if ion_num == int(pdos_matrix[j][0]):
                        ion_pdos_data_row = [band_num, band_en, band_occ, spin_component]
                        for i in range(1,lm_orb_tot+1):
                            ion_pdos_data_row.append(float(pdos_matrix[j][i])/norm_const)
                ion_pdos_data.append(ion_pdos_data_row)

    outdata = '#band      energy (eV)  occ spin'
    for orbital in lm_orbitals:
        outdata += '    %5s'%orbital
    outdata += '\n'
    for row in ion_pdos_data:
        outdata += '%5d    % 13.8f    %1d    %1d'%(row[0], row[1], row[2], row[3])
        for i in range(4, len(row)):
            outdata += '    %1.3f'%row[i]
        outdata += '\n'
    
    if plot_data:
        import matplotlib.pyplot as plt
        import numpy as np
        ion_pdos_plot = []
        for row in ion_pdos_data:
            ion_pdos_plot_row = []
            ion_pdos_plot_row.append(row[1]) # energy
            ion_pdos_plot_row.append(row[2]) # occ
            ion_pdos_plot_row.append(row[3]) # spin
            ion_pdos_plot_row.append(row[4]) # l = 0 state
            ion_pdos_plot_row.append(row[2]*row[4]) # l = 0 occ state
            ion_pdos_plot_row.append(row[5]+row[6]+row[7]) # l = 1 state
            ion_pdos_plot_row.append(row[2]*(row[5]+row[6]+row[7])) # l = 1 occ state
            ion_pdos_plot.append(ion_pdos_plot_row)
        ion_pdos_plot = np.array(ion_pdos_plot)

        def plot_gauss_smearing(x_energies, y_dos, sigma):
            x_min = min(x_energies)
            x_max = max(x_energies)
            x_plot = np.linspace(x_min, x_max, 10000)
            y_plot = np.zeros(10000)
            for i in range(len(x_energies)):
                mu = x_energies[i]
                y_amp = y_dos[i]
                y_plot = y_plot + y_amp*gauss(x_plot, mu, sigma)
            return x_plot, y_plot

        band_energies, ion_pdos_s_up = plot_gauss_smearing(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,3], smearing)
        band_energies, ion_pdos_s_up_occ = plot_gauss_smearing(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,4], smearing)
        plt.plot(band_energies, ion_pdos_s_up, label='s', color='red')
        plt.fill_between(band_energies, ion_pdos_s_up_occ, color='red')

        band_energies, ion_pdos_s_dwn = plot_gauss_smearing(ion_pdos_plot[bands_total:,0], ion_pdos_plot[bands_total:,3], smearing)
        band_energies, ion_pdos_s_dwn_occ = plot_gauss_smearing(ion_pdos_plot[bands_total:,0], ion_pdos_plot[bands_total:,4], smearing)
        plt.plot(band_energies, -ion_pdos_s_dwn, color='red')
        plt.fill_between(band_energies, -ion_pdos_s_dwn_occ, color='red')

        band_energies, ion_pdos_p_up = plot_gauss_smearing(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,5], smearing)
        band_energies, ion_pdos_p_up_occ = plot_gauss_smearing(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,6], smearing)
        plt.plot(band_energies, ion_pdos_p_up, label='p', color='green')
        plt.fill_between(band_energies, ion_pdos_p_up_occ, color='green')

        band_energies, ion_pdos_p_dwn = plot_gauss_smearing(ion_pdos_plot[bands_total:,0], ion_pdos_plot[bands_total:,5], smearing)
        band_energies, ion_pdos_p_dwn_occ = plot_gauss_smearing(ion_pdos_plot[bands_total:,0], ion_pdos_plot[bands_total:,6], smearing)
        plt.plot(band_energies, -ion_pdos_p_dwn, color='green')
        plt.fill_between(band_energies, -ion_pdos_p_dwn_occ, color='green')
        
        #plt.plot(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,3], label='s', color='red', linestyle='none', marker='o', markersize=1)
        #plt.vlines(ion_pdos_plot[:bands_total,0], 0*ion_pdos_plot[:bands_total,4], ion_pdos_plot[:bands_total,4], color='red')
        #plt.fill_between(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,3], where=ion_pdos_plot[:bands_total,1]>0, color='red')
        
        #plt.plot(ion_pdos_plot[bands_total:,0], -ion_pdos_plot[bands_total:,3], label='s', color='red', linestyle='none', markersize=1)
        #plt.vlines(ion_pdos_plot[bands_total:,0], 0*ion_pdos_plot[bands_total:,4], -ion_pdos_plot[bands_total:,4], color='red')
        #plt.fill_between(ion_pdos_plot[bands_total:,0], -ion_pdos_plot[bands_total:,3], where=ion_pdos_plot[bands_total:,1]>0, color='red')
        
        #plt.plot(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,5], label='p', color='green', linestyle='none', markersize=1)
        #plt.vlines(ion_pdos_plot[:bands_total,0], 0*ion_pdos_plot[:bands_total,6], ion_pdos_plot[:bands_total,6], color='green')
        #plt.fill_between(ion_pdos_plot[:bands_total,0], ion_pdos_plot[:bands_total,4], where=ion_pdos_plot[:bands_total,1]>0, color='green')

        #plt.plot(ion_pdos_plot[bands_total:,0], -ion_pdos_plot[bands_total:,5], label='p', color='green', linestyle='none', markersize=1)
        #plt.vlines(ion_pdos_plot[bands_total:,0], 0*ion_pdos_plot[bands_total:,5], -ion_pdos_plot[bands_total:,5], color='green')
        #plt.fill_between(ion_pdos_plot[bands_total:,0], -ion_pdos_plot[bands_total:,4], where=ion_pdos_plot[bands_total:,1]>0, color='green')

        #plt.xlim([-5,4])
        plt.title("PDOS of ion no. %d"%ion_num)
        plt.xlabel("Energy (eV)")
        plt.yticks([])
        plt.legend()
        plt.savefig(file_out_name + '.png')
        plt.show()
    # Writing file
    outfile = open(file_out_name,'w')
    outfile.write(outdata)
    outfile.close()
    return "Successfull termination."

def gauss(x, mu, sigma):
    import numpy as np
    var = sigma**2
    return (1/np.sqrt(2*np.pi*var))*np.exp(-(x-mu)**2/(2*var))

def bohr2A(val):
    return 0.529177249*val

#Matrix determinant
def mat_det(m):
    n=len(m)
    if n != len(m[0]):
        raise ValueError("Determinant called for a non-square matrix.")
    if n == 1:
        return m[0][0]
    if n == 2:
        det2x2=(m[0][0]*m[1][1])-(m[0][1]*m[1][0])
        return det2x2
    else:
        detnxn = 0
        for i in range(n):
            column = []
            for y in range(1, n):
                row = []
                for x in range(i):
                    row.append(m[y][x])
                for x in range(i+1,n):
                    row.append(m[y][x])
                column.append(row)
            if i % 2 == 0:
                detnxn+=m[0][i]*mat_det(column)
            else:
                detnxn-=m[0][i]*mat_det(column)
    return detnxn

#Matrix Transpose
def mat_tran(m):
    leny = len(m)
    lenx = len(m[0])
    mT = []
    for x in range(lenx):
        row = []
        for y in range(leny):
            row.append(m[y][x])
        mT.append(row)
    return mT

#Matrix Adjoint
def mat_adj(m):
    n=len(m)
    if n != len(m[0]):
        raise ValueError("Adjoint called for a non-square matrix.")
    adjM = []
    for j in range(n):
        adjrow = []
        for i in range(n):
            minor = []
            for y in range(j):
                row = []
                for x in range(i):
                    row.append(m[y][x])
                for x in range(i+1,n):
                    row.append(m[y][x])
                minor.append(row)
            for y in range(j+1,n):
                row = []
                for x in range(i):
                    row.append(m[y][x])
                for x in range(i+1,n):
                    row.append(m[y][x])
                minor.append(row)
            cofactor = (-1)**(i+j)
            adjrow.append(cofactor*mat_det(minor))
        adjM.append(adjrow)
    return mat_tran(adjM)

def mat_inv(m):
    n=len(m)
    if n != len(m[0]):
        raise ValueError("Inverse called for a non-square matrix.")
    Minv = []
    Madj = mat_adj(m)
    Mdet = mat_det(m)
    for j in range(len(m)):
        Mrowj = []
        for i in range(len(m)):
            Mrowj.append(Madj[j][i]/Mdet)
        Minv.append(Mrowj)
    return Minv

if __name__ == "__main__":
    import sys
    sys.exit(main())