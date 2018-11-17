#!/usr/bin/env python

import argparse

from kaldi.util.table import SequentialLatticeReader, LatticeWriter
from kaldi.fstext import LatticeVectorFst, LatticeWeight, LatticeArc
from kaldi.fstext import LatticeVectorFstArcIterator, LatticeVectorFstStateIterator
from kaldi.fstext import compose, properties

def RemoveCTCBlankFromLattice(lat, blank):
    # Get mapping from output symbols from the output labels in the input lattice to states
    symbol2state = {}
    for s in lat.states():
        for a in lat.arcs(s):
            o = a.olabel;
            if o != blank and o not in symbol2state:
                symbol2state[o] = len(symbol2state)+1

    # Create composition lattice, such that output = compose(input, C)
    C = LatticeVectorFst()
    for n in range(len(symbol2state)+1):
        C.set_final(C.add_state(), LatticeWeight.one())
    C.set_start(0)

    # Self-loop in the blank state
    C.add_arc(0, LatticeArc(blank, 0, LatticeWeight.one(), 0))

    for label, state_id in symbol2state.items():
        # Arc from the initial state to the symbol's state
        C.add_arc(0, LatticeArc(label, label, LatticeWeight.one(), state_id))
        # Self-loop in the symbol's state
        C.add_arc(state_id, LatticeArc(label, 0, LatticeWeight.one(), state_id))
        # Arc back to the blank state
        C.add_arc(state_id, LatticeArc(blank, 0, LatticeWeight.one(), 0))
        # Arc to all other symbol states
        for label2, state_id2 in symbol2state.items():
            if label != label2:
                C.add_arc(state_id, LatticeArc(label2, label2, LatticeWeight.one(), state_id2))

    # Compute output lattice
    return compose(lat, C)


### Main block called when run from command line ###
if __name__ == '__main__':

    ### Parse input arguments ###
    def raise_(ex):
        """Raise that works within lambda functions."""
        raise ex
    parser = argparse.ArgumentParser(description='Remove CTC blank symbols from the output labels of Kaldi lattices.')
    parser.add_argument('blank', type=lambda x: int(x) if int(x) > 0 else raise_(ValueError),
                        help='Blank symbol index.')
    parser.add_argument('rspecifier', type=str,
                        help='Input lattice sequence.')
    parser.add_argument('wspecifier', type=str,
                        help='Output lattice sequence.')
    args = parser.parse_args()

    ### Initialize output writer ###
    outlat = LatticeWriter(args.wspecifier)

    ### Loop through input lattices ###
    for lat_key, lat in SequentialLatticeReader(args.rspecifier):

        ### Make sure that lattice complies with all assumptions ##
        prop = lat.properties(properties.ACCEPTOR | properties.ACYCLIC, True)
        if prop & properties.ACCEPTOR != properties.ACCEPTOR:
            raise Exception('Lattice '+lat_key+' is not an acceptor.')
        if prop & properties.ACYCLIC != properties.ACYCLIC:
            raise Exception('Lattice '+lat_key+' is not an acyclic.')

        ### Remove blank ###
        out = RemoveCTCBlankFromLattice(lat, args.blank)

        ### Write to output ###
        outlat.write(lat_key, out)
