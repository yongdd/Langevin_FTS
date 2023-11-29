/*----------------------------------------------------------
* This class defines class for propagators
*-----------------------------------------------------------*/

#ifndef PROPAGATORS_H_
#define PROPAGATORS_H_

#include <string>
#include <vector>
#include <map>
#include "Polymer.h"

struct ComputationEdge{
    int max_n_segment;                                    // the maximum segment number
    std::string monomer_type;                             // monomer_type
    std::vector<std::tuple<std::string, int, int>> deps;  // tuple <key, n_segment, n_repeated>
    int height;                                           // height of propagator (height of tree data Structure)
};
struct ComputationBlock{
    std::string monomer_type;  // monomer_type

    // When the 'aggregate_propagator_computation' is on, original block can be sliced to smaller block pieces.
    // For example, suppose one block is composed of 5'A' monomers. -> -A-A-A-A-A-
    // And this block is sliced into 3'A' monomers and 2'A' monomers -> -A-A-A-,  -A-A-
    // For the first slice, n_segment_allocated, n_segment_offset, and n_segment_original are 3, 0, and 5, respectively.
    // For the second slice, n_segment_allocated, n_segment_offset, and n_segment_original are 2, 3, and 5, respectively.
    // If the 'aggregate_propagator_computation' is off, original block is not sliced to smaller block pieces.
    // In this case, n_segment_allocated, n_segment_offset, and n_segment_original are 5, 0, and 5, respectively.

    int n_segment_allocated;
    int n_segment_offset;      
    int n_segment_original;
    std::vector<std::tuple<int ,int>> v_u; // node pair <polymer id, v, u>
};

/* This stucture defines comparison function for branched key */
struct ComparePropagatorKey
{
    bool operator()(const std::string& str1, const std::string& str2);
};

// class Propagators
// {
// private:
//     std::string model_name; // "continuous": continuous standard Gaussian model
//                             // "discrete": discrete bead-spring model
                            
//     bool aggregate_propagator_computation; // compute multiple propagators using property of linearity of the diffusion equation.

//     // set{key: (polymer id, dep_v, dep_u) (assert(dep_v <= dep_u))}
//     std::map<std::tuple<int, std::string, std::string>, ComputationBlock> essential_blocks;

//     // dictionary{key:non-duplicated unique propagator_codes, value: ComputationEdge}
//     std::map<std::string, ComputationEdge, ComparePropagatorKey> essential_propagator_codes; 

//     // Add new key. if it already exists and 'new_n_segment' is larger than 'max_n_segment', update it.
//     void update_essential_propagator_code(std::map<std::string, ComputationEdge, ComparePropagatorKey>& essential_propagator_codes, std::string new_key, int new_n_segment);

//     // Superpose propagators
//     std::map<std::string, ComputationBlock> superpose_propagator_common             (std::map<std::string, ComputationBlock> remaining_keys, int minimum_n_segment);
//     std::map<std::string, ComputationBlock> superpose_propagator_of_continuous_chain(std::map<std::string, ComputationBlock> u_map);
//     std::map<std::string, ComputationBlock> superpose_propagator_of_discrete_chain  (std::map<std::string, ComputationBlock> u_map);

// public:

//     Propagators();
//     ~Propagators() {};

//     // Get information of essential propagators and blocks
//     bool is_using_propagator_aggregation() const;
//     int get_n_essential_propagator_codes() const;
//     std::map<std::string, ComputationEdge, ComparePropagatorKey>& get_essential_propagator_codes(); 
//     ComputationEdge& get_essential_propagator_code(std::string key);
//     std::map<std::tuple<int, std::string, std::string>, ComputationBlock>& get_essential_blocks(); 
//     ComputationBlock& get_essential_block(std::tuple<int, std::string, std::string> key);

//     // Display
//     void display_propagators() const;
//     void display_blocks() const;
//     void display_sub_propagators() const;
// };
#endif
