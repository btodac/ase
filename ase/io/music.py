""" Music mol file writer
"""


from ase.atoms import Atoms

# Element dict
symbol2name = {'Ac':'Actinium',
    'Ag':'Silver',
    'Al':'Aluminum',
    'Am':'Americium',
    'Ar':'Argon',
    'As':'Arsenic',
    'At':'Astatine',
    'Au':'Gold',
    'Ba':'Barium',
    'Bo':'Boron',
    'Be':'Beryllium',
    'Bh':'Bohrium',
    'Bi':'Bismuth',
    'Bk':'Berkelium',
    'Br':'Bromine',
    'Ca':'Calcium',
    'Cd':'Cadmium',
    'C' :'Carbon',
    'Ce':'Cerium',
    'Cf':'Californium',
    'Cl':'Chlorine',
    'Cn':'Copernicium',
    'Cm':'Curium',
    'Co':'Cobalt',
    'Cr':'Chromium',
    'Cs':'Cesium',
    'Cu':'Copper',
    'Ds':'Darmstadtium',
    'Db':'Dubnium',
    'Dy':'Dysprosium',
    'Er':'Erbium',
    'Es':'Einsteinium',
    'Eu':'Europium',
    'Fm':'Fermium',
    'Fl':'Flerovium',
    'F':'Fluorine',
    'Fe':'Iron',
    'Fr':'Francium',
    'Ga':'Gallium',
    'Gd':'Gadolinium',
    'Ge':'Germanium',
    'H':'Hydrogen',
    'He':'Helium',
    'Hg':'Mercury',
    'Hf':'Hafnium',
    'Ho':'Holmium',
    'Hs':'Hassium',
    'I' :'Iodine',
    'In':'Indium',
    'Ir':'Iridium',
    'K':'Potassium',
    'Kr':'Krypton',
    'La':'Lanthanum',
    'Li':'Lithium',
    'Lr':'Lawrencium',
    'Lv':'Livermorium',
    'Lu':'Lutetium',
    'Md':'Mendelevium',
    'Mg':'Magnesium',
    'Mn':'Manganese',
    'Mt':'Meitnerium',
    'Mo':'Molybdenum',
    'Mc':'Moscovium',
    'N' :'Nitrogen',
    'Na':'Sodium',
    'Nb':'Niobium',
    'Nd':'Neodymium',
    'Ne':'Neon',
    'Ni':'Nickel',
    'Nh':'Nihonium',
    'No':'Nobelium',
    'Np':'Neptunium',
    'O' :'Oxygen',
    'Og':'Oganessonm',
    'Os':'Osmium',
    'Pb':'Lead',
    'P' :'Phosphorus',
    'Pa':'Protactinium',
    'Pd':'Palladium',
    'Po':'Polonium',
    'Pr':'Praseodymium',
    'Pm':'Promethium',
    'Pt':'Platinum',
    'Pu':'Plutonium',
    'Ra':'Radium',
    'Rb':'Rubidium',
    'Re':'Rhenium',
    'Rf':'Rutherfordium',
    'Rg':'Roentgenium',
    'Rh':'Rhodium',
    'Rn':'Radon',
    'Ru':'Ruthenium',
    'S' :'Sulfur',
    'Sb':'Antimony',
    'Sc':'Scandium',
    'Sm':'Samarium',
    'Sg':'Seaborgium',
    'Se':'Selenium',
    'Si':'Silicon',
    'Sn':'Tin',
    'Sr':'Strontium',
    'Ta':'Tantalum',
    'Tb':'Terbium',
    'Tc':'Technetium',
    'Te':'Tellurium',
    'Th':'Thorium',
    'Ts':'Tennessine',
    'Tl':'Thallium',
    'Ti':'Titanium',
    'Tm':'Thulium',
    'W' :'Tungsten',
    'U' :'Uranium',
    'V' :'Vanadium',
    'Xe':'Xenon',
    'Y' :'Yttrium',
    'Yb':'Ytterbium',
    'Zn':'Zinc',
    'Zr':'Zirconium'}

# Inverse dict
name2symbol = {v: k for k, v in symbol2name.items()}


def read_music(fileobj):
    """ Method to read Music mol files into Atoms object
    """

    lines = fileobj.readlines()
    L1 = lines[0].split()
    mol_name = L1[1]
    del(lines[:3])
    L1 = lines[0].split()
    del(lines[0])
    natoms = int(L1[0])
    positions = []
    symbols = []
    charges = []
    for line in lines[:natoms]:
            N, x, y, z, element, Q, s, t = line.split()[:8]
            symbol = name2symbol[element]
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
            charges.append(float(Q))
    return Atoms(symbols=symbols, positions=positions,charges=charges)

def write_music(fileobj, images, name='None'):
    """ Method to  write Music mol files
    """
    # Get first atoms object if images is a list
    if isinstance(images, (list,)):
        atoms = images[0]
    else:
        atoms = images
    # extract the data
    symbols = atoms.get_chemical_symbols()
    natoms = len(symbols)
    coords = atoms.get_positions()
    charges = atoms.get_initial_charges()
    cell = atoms.get_cell_lengths_and_angles()
    celldisp = atoms.get_celldisp()
    # Write the header
    fileobj.write(' Molecule_name: %s\n' % name)
    fileobj.write('\n')
    fileobj.write('  Coord_Info: Listed Cartesian None\n')
    fileobj.write('         %d\n' % natoms)
    # write the atoms information
    for i,s in enumerate(symbols):
        x, y, z = coords[i][:]
        Q = charges[i]
        s = symbol2name[s]
        fileobj.write('%5s %9.5f %9.5f %9.5f   %-13s %6.2f  0  0\n' 
                % (i+1, x, y, z, s, Q))
        # Num, x, y ,z coords, atom name (element name!!!), charge
    fileobj.write('\n\n\n')
    fileobj.write('  Fundcell_Info: Listed\n')
    fileobj.write('  %13.5f  %13.5f  %13.5f\n' % (cell[0], cell[1], cell[2]))
    fileobj.write('  %13.5f  %13.5f  %13.5f\n' % (cell[3], cell[4], cell[5]))
    fileobj.write('  %13.5f  %13.5f  %13.5f\n' 
                    % (celldisp[0], celldisp[1], celldisp[2]))
    fileobj.write('  %13.5f  %13.5f  %13.5f\n' % (cell[0], cell[1], cell[2]))




