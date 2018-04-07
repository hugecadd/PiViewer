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

import openbabel as ob
import pybel
import numpy as np

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
def find_PiPi(pdb_file, lig_name, centroid_distance=5.0, dih_parallel=25, dih_tshape=80, verbose=1):
    """
    Find Pi-Pi interactions around the specified ligand residue from the pdb file.
    :param pdb_file: path of the target file in PDB format.
    :param lig_name: ligand residue name.
    :param centroid_distance: Max ring centroid distance
    :param dih_parallel: Max dihedral (parallel)
    :param dih_tshape: Min dihedral (T-shaped)
    :return: number of Pi-Pi interactions found
    """
    # Get ligand residue and print its name.
    ligAtomList = []
    ligAtomIdList = []
    mol = next(pybel.readfile('pdb', pdb_file))
    if verbose: print("A total of %s residues" % mol.OBMol.NumResidues())
    lig = None
    for res in ob.OBResidueIter(mol.OBMol):
        # print res.GetName()
        if res.GetName() == lig_name:
            lig = res
            if verbose: print("Ligand residue name is:", lig.GetName())
            break
    if not lig:
        if verbose: print("No ligand residue %s found, please confirm." % lig_name)
        return -1
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
                    if verbose: print("ligand ring_ID: ", ring.ring_id, end=' ')
                    if ring.IsAromatic():
                        if verbose: print("aromatic")
                        ligAroRingList.append(ring)
                    else:
                        if verbose: print("saturated")
    for ring in mol.sssr:
        if ring.ring_id not in ligRingIdList:
            recRingList.append(ring)
            if ring.IsAromatic():
                recAroRingList.append(ring)
    if verbose: print("\nReceptor has ", len(recRingList), " rings,", end=' ')
    if verbose: print(" has ", len(recAroRingList), " aromatic rings.")

    # Find and show the rings
    ligRingCenter = ob.vector3()
    recRingCenter = ob.vector3()
    ligNorm1 = ob.vector3()
    ligNorm2 = ob.vector3()
    recNorm1 = ob.vector3()
    recNorm2 = ob.vector3()
    count = 0
    lig_ring_index = 0
    for ligRing in ligAroRingList:
        lig_ring_index += 1
        ligRing.findCenterAndNormal(ligRingCenter, ligNorm1, ligNorm2)
        rec_ring_index = 0
        for recRing in recAroRingList:
            rec_ring_index += 1
            recRing.findCenterAndNormal(recRingCenter, recNorm1, recNorm2)
            dist = ligRingCenter.distSq(recRingCenter) ** 0.5
            angle = vecAngle(ligNorm1, recNorm1)
            if (dist < centroid_distance and (angle < dih_parallel or angle > dih_tshape)):  # the criteria
                count += 1
                if verbose: print("Pi-Pi ring pairs: %3s,%3s  Angle(deg.): %5.2f  Distance(A): %.2f" % (recRing.ring_id, ligRing.ring_id, angle, dist))
    if verbose: print("Total Pi-Pi interactions:", count)
    return count


if __name__ == '__main__':
    pdb_file = r'C:\CloudStation\Epicat\Git\PiViewer\1ACJ.pdb'
    lig_name = 'THA'
    find_PiPi(pdb_file, lig_name, 5.0, 25, 80)
