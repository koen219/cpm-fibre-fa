from ecmgen.networks import random_network, single_strand, single_spring, laminin
from ecmgen.network import Network

import unittest


class TestAddingNetwork(unittest.TestCase):
    def test_empty_networks(self):
        net1 = single_strand(200, 200, 100, 100, 0, 9, 50)
        net2 = single_strand(200, 200, 100, 100, 3.1415 * 0.5, 9, 50)

        net3 = net1 + net2
        self.assertIsInstance(net3, Network)
        

    def test_adding_two_strands(self):
        net1 = single_strand(200, 200, 100, 100, 0, 9, 50)
        net2 = single_strand(200, 200, 100, 100, 3.1415 * 0.5, 9, 50)

        net3 = net1 + net2

        self.assertEqual(
            len(net3.beads_positions),
            len(net1.beads_positions) + len(net2.beads_positions),
        )
        self.assertEqual(
            len(net3.bonds_groups), len(net1.bonds_groups) + len(net2.bonds_groups)
        )
        self.assertEqual(
            len(net3.angle_groups), len(net1.angle_groups) + len(net2.angle_groups)
        )

import itertools
class TestLaminin(unittest.TestCase):
    def test_addingLaminin(self):
        network = random_network(
            sizex=200,
            sizey=200,
            number_of_beads_per_strand=9,
            number_of_strands=10,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=0,
            crosslink_bin_size=1 / 3,
            seed=10,
        )
        laminin(sizex=200,sizey=200,amount_of_laminin=10,network=network,seed=1)
        self.assertSetEqual(
            set(network.bonds_types),
            set(["polymer", "laminin"])
        )  
        self.assertSetEqual(
            set(network.beads_types),
            set(["free", "boundary"])
        )  
        laminin_ids = [k for (k,typ) in enumerate(network.bonds_types) if typ == 'laminin']
        self.assertEqual(len(laminin_ids), 10)

        self.assertEqual(
            len(network.beads_positions),
            len(network.beads_types)
        )

        self.assertEqual(
            len(network.bonds_groups),
            len(network.bonds_types)
        )

    def test_particleids(self):
        network = random_network(
            sizex=200,
            sizey=200,
            number_of_beads_per_strand=9,
            number_of_strands=10,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=0,
            crosslink_bin_size=1 / 3,
            seed=10,
        )
        laminin(sizex=200,sizey=200,amount_of_laminin=10,network=network,seed=1)

        self.assertLess(
            len([bond[0] for bond in network.bonds_groups]),
            len(network.beads_positions),
        )
        self.assertLess(
            len([bond[1] for bond in network.bonds_groups]),
            len(network.beads_positions),
        )
        
        self.assertLess(
            max(list(itertools.chain(*network.bonds_groups))),
            len(network.beads_positions)
        )


class TestStringMethods(unittest.TestCase):
    def test_creation(self):
        network = random_network(
            sizex=200,
            sizey=200,
            number_of_beads_per_strand=9,
            number_of_strands=100,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=50,
            crosslink_bin_size=1 / 3,
            seed=10,
        )

    def test_crosslinks_creation(self):
        network = random_network(
            sizex=200,
            sizey=200,
            number_of_beads_per_strand=9,
            number_of_strands=100,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=50,
            crosslink_bin_size=1 / 3,
            seed=10,
        )
        for bondtype in network.bonds_types:
            if bondtype == "polymer":
                continue
            self.assertLessEqual(network.details_of_bondtypes[bondtype]["r0"], 1.0)

    def test_singleStrand(self):
        beads = 9
        network = single_strand(
            200, 200, 50, 50, 3.141592653589793, beads, 1 * (beads - 1), None
        )
        for i in range(beads):
            self.assertEqual(network.beads_positions[i][0], 50 + i)
            self.assertEqual(network.beads_positions[i][1], 50)

    def test_singleSpring2(self):
        network = single_spring(
            sizex=200,
            sizey=200,
            number_of_beads_per_strand=3,
            number_of_strands=11,
            contour_length_of_strand=6.25,
            seed=None,
        )

    #    def test_crosslinks_creation(self):
    #        network = random_network(
    #            sizex=200,
    #            sizey=200,
    #            number_of_beads_per_strand=9,
    #            number_of_strands=100,
    #            contour_length_of_strand=50,
    #            crosslink_max_r=1.0,
    #            maximal_number_of_initial_crosslinks=50,
    #            crosslink_bin_size=1 / 3,
    #            seed=10,
    #        )
    #        for bondtype in network.bonds_types:
    #            if bondtype == "polymer":
    #                continue
    #            self.assertLessEqual(network.details_of_bondtypes[bondtype]['r0'], 1.0)
    #        self.assertEqual(sum([1 for x in network.bonds_types if x == "polymer"]), 11 * 2)
    #
    #        self.assertEqual(len(network.bonds_groups), 32)
    #        self.assertEqual(list(network.bonds_groups[0]), [0,1])
    #        self.assertEqual(list(network.bonds_groups[1]), [1,2])
    #        self.assertEqual(list(network.bonds_groups[2]), [3,4])
    #        self.assertEqual(list(network.bonds_groups[3]), [4,5])
    #        self.assertEqual(list(network.bonds_groups[4]), [6,7])
    #        self.assertEqual(list(network.bonds_groups[5]), [7,8])
    #        self.assertEqual(list(network.bonds_groups[6]), [9,10])
    #        self.assertEqual(list(network.bonds_groups[7]), [10,11])
    #        self.assertEqual(list(network.bonds_groups[2*11]), [2,3])
    #        self.assertEqual(list(network.bonds_groups[2*11+1]), [5,6])
    #        self.assertEqual(list(network.bonds_groups[2*11+2]), [8,9])
    #
    #        for bond in network.bonds_groups:
    #            self.assertLess(bond[0] , len(network.beads_positions))
    #            self.assertLess(bond[1] , len(network.beads_positions))

    def test_boundary_particles(self):
        network = random_network(
            sizex=50,
            sizey=50,
            number_of_beads_per_strand=9,
            number_of_strands=100,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=50,
            crosslink_bin_size=1 / 3,
            seed=10,
            fix_boundary=True,
        )
        self.assertIn("boundary", network.beads_types)

        network = random_network(
            sizex=50,
            sizey=50,
            number_of_beads_per_strand=9,
            number_of_strands=100,
            contour_length_of_strand=50,
            crosslink_max_r=1.0,
            maximal_number_of_initial_crosslinks=50,
            crosslink_bin_size=1 / 3,
            seed=10,
            fix_boundary=False,
        )
        self.assertNotIn("boundary", network.beads_types)


if __name__ == "__main__":
    unittest.main()
