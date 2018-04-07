# -*- coding: utf-8 -*-
"""
PiViewer: an open-source tool for automated detection and display of pi-pi interactions.
Created by Hu Ge <gehuchina@gmail.com>
"""

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# #
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

from pymol import cmd
from Tkinter import *
import tkFileDialog  # not included in TKinter standard
import openbabel as ob
import pybel
import numpy as np
import os


# Define vecAngle for the degree between two vectors.
def vecAngle(vec1, vec2):
    '''
    vec1,vec2 are ob.vector3 objects,
    return the angle between them in degree.
    '''
    dotprod = vec1.GetX() * vec2.GetX() + vec1.GetY() * vec2.GetY() + vec1.GetZ() * vec2.GetZ()
    deg = np.arccos(dotprod) * 180 / np.pi
    if deg > 90:
        deg = 180 - deg
    return deg


# The main PiPi viewer function
def PiPi(pdb_file, lig_name, centroid_distance, dih_parallel, dih_tshape):
    # Get ligand residue and print its name.
    ligAtomList = []
    ligAtomIdList = []
    mol = pybel.readfile('pdb', pdb_file).next()
    print "A total of %s residues" % mol.OBMol.NumResidues()
    lig = None
    for res in ob.OBResidueIter(mol.OBMol):
        # print res.GetName()
        if res.GetName() == lig_name:
            lig = res
            print "Ligand residue name is: ", lig.GetName()
            break
    if not lig:
        print "No ligand residue %s found, please confirm." % lig_name
        return 0
    else:
        for atom in ob.OBResidueAtomIter(lig):
            # print atom.GetIdx()
            ligAtomList.append(atom)
            ligAtomIdList.append(atom.GetIdx())

    # Set ring_id
    i = 0
    for ring in mol.sssr:
        ring.ring_id = i
        i += 1
        # print ring.ring_id

    # Determine which rings are from ligand.
    ligRingList = []
    ligAroRingList = []
    ligRingIdList = []
    recRingList = []
    recAroRingList = []
    for ring in mol.sssr:
        for atom in ligAtomList:
            if ring.IsMember(atom):
                if ring not in ligRingList:
                    ligRingList.append(ring)
                    ligRingIdList.append(ring.ring_id)
                    print "ligand ring_ID: ", ring.ring_id,
                    if ring.IsAromatic():
                        print "Aromatic"
                        ligAroRingList.append(ring)
                    else:
                        print "Saturated"
    for ring in mol.sssr:
        if ring.ring_id not in ligRingIdList:
            recRingList.append(ring)
            if ring.IsAromatic():
                recAroRingList.append(ring)
    print "\nReceptor has ", len(recRingList), " rings,",
    print " has ", len(recAroRingList), " aromatic rings."

    # Find and show the rings
    ligRingCenter = ob.vector3()
    recRingCenter = ob.vector3()
    ligNorm1 = ob.vector3()
    ligNorm2 = ob.vector3()
    recNorm1 = ob.vector3()
    recNorm2 = ob.vector3()
    coord1 = []
    coord2 = []
    pair = []  # Store the names of the objects
    pairList = []
    psObjectList = []
    i = 0
    for ligRing in ligAroRingList:
        ligRing.findCenterAndNormal(ligRingCenter, ligNorm1, ligNorm2)
        coord1 = [ligRingCenter.GetX(), ligRingCenter.GetY(), ligRingCenter.GetZ()]
        # Create a pseudoatom for the centroid of each ring in the ligand
        objectName1 = 'psCenter%.2d' % i
        i += 1
        cmd.pseudoatom(object=objectName1, pos=coord1)
        psObjectList.append(objectName1)
        for recRing in recAroRingList:
            recRing.findCenterAndNormal(recRingCenter, recNorm1, recNorm2)
            dist = ligRingCenter.distSq(recRingCenter)
            angle = vecAngle(ligNorm1, recNorm1)
            if (dist ** 0.5 < centroid_distance and (angle < dih_parallel or angle > dih_tshape)):  # the criteria
                coord2 = [recRingCenter.GetX(), recRingCenter.GetY(), recRingCenter.GetZ()]
                # Create pseudoatom for each lig ring center
                objectName2 = 'psCenter%.2d' % i
                i += 1
                cmd.pseudoatom(object=objectName2, pos=coord2)
                psObjectList.append(objectName2)
                # pair=[coord1,coord2]
                pair = [objectName1, objectName2]
                # print "%.4f,%.4f,%.4f" % (pair[1][0],pair[1][1],pair[1][2])
                pairList.append(pair)
                cmd.distance('dist-' + objectName1 + '-' + objectName2, pair[0], pair[1])
                print objectName1, objectName2, " angle is : %.2f" % angle


# Define the dialog for parameters
def setParamDialog(app):
    top = Tk()
    top.wm_title('PiViewer v1.2 for PyMOL')
    # Labels
    L1 = Label(top, text="Ligand residue name")
    L2 = Label(top, text="Max ring centroid distance")
    L3 = Label(top, text="Max dihedral (parallel)")
    L4 = Label(top, text="Min dihedral (T-shaped)")
    L1.grid(row=0, column=0, sticky=W)
    L2.grid(row=1, column=0, sticky=W)
    L3.grid(row=2, column=0, sticky=W)
    L4.grid(row=3, column=0, sticky=W)
    # Text entries and default values
    E1 = Entry(top, bd=5)
    E2 = Entry(top, bd=5)
    E3 = Entry(top, bd=5)
    E4 = Entry(top, bd=5)
    # Modify the Default values here
    E1.insert(0, "LIG")
    E2.insert(0, "5.0")
    E3.insert(0, "25.0")
    E4.insert(0, "80.0")
    E1.grid(row=0, column=1)
    E2.grid(row=1, column=1)
    E3.grid(row=2, column=1)
    E4.grid(row=3, column=1)

    # When clicked, get the input and launch the main Pi-Pi viewer function - PiPi
    def B1Click():
        lig_name = E1.get()
        centroid_distance = float(E2.get())
        dih_parallel = float(E3.get())
        dih_tshape = float(E4.get())
        pdb_file = tkFileDialog.askopenfilename(
            filetypes=[("PDB(Protein Data Bank files", "*.pdb"), ("PDB(Protein Data Bank files", "*.ent"),
                       ('All files', '*')], title='Please select the target PDB file')
        pdb_file = str(os.path.normpath(pdb_file))
        print lig_name, centroid_distance, dih_parallel, dih_tshape
        print pdb_file
        PiPi(pdb_file, lig_name, centroid_distance, dih_parallel, dih_tshape)
        top.destroy()

    B1 = Button(top, bd=2, pady=5, text="Run", command=B1Click)
    B2 = Button(top, bd=2, pady=5, text="Cancel", command=top.destroy)
    B1.grid(row=4, column=0, sticky=W + E)
    B2.grid(row=4, column=1, sticky=W + E)
    top.mainloop()


def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                             'PiViewer',
                             label='PiViewer',
                             command=lambda s=self: setParamDialog(s))

cmd.extend('PiPi', PiPi)
