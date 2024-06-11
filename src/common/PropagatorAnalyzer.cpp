#include <iostream>
#include <cctype>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stack>
#include <set>

#include "PropagatorAnalyzer.h"
#include "Molecules.h"
#include "Polymer.h"
#include "Exception.h"

bool ComparePropagatorKey::operator()(const std::string& str1, const std::string& str2)
{
    // First compare heights
    int height_str1 = PropagatorCode::get_height_from_key(str1);
    int height_str2 = PropagatorCode::get_height_from_key(str2);

    if (height_str1 < height_str2)
        return true;
    else if(height_str1 > height_str2)
        return false;

    // Second compare their strings
    return str1 > str2;
}

PropagatorAnalyzer::PropagatorAnalyzer(Molecules* molecules, bool aggregate_propagator_computation)
{
    if(molecules->get_n_polymer_types() == 0)
        throw_with_line_number("There is no chain. Add polymers first.");

    this->aggregate_propagator_computation = aggregate_propagator_computation;
    this->model_name = molecules->get_model_name();
    for(int p=0; p<molecules->get_n_polymer_types();p++)
    {
        add_polymer(molecules->get_polymer(p), p);
    }
}
void PropagatorAnalyzer::add_polymer(Polymer& pc, int polymer_id)
{
    // Temporary map for the new polymer
    std::map<std::string, std::map<std::string, ComputationBlock>> computation_blocks_new_polymer;
    std::map<std::tuple<int, int>, std::string> v_u_to_right_key;

    // Find computation_blocks in new_polymer
    std::vector<Block> blocks = pc.get_blocks();
    for(size_t b=0; b<blocks.size(); b++)
    {
        int v = blocks[b].v;
        int u = blocks[b].u;
        std::string dep_v = pc.get_propagator_key(v, u);
        std::string dep_u = pc.get_propagator_key(u, v);

        if (dep_v < dep_u){
            dep_v.swap(dep_u);
            std::swap(v,u);
        }

        computation_blocks_new_polymer[dep_v][dep_u].monomer_type = blocks[b].monomer_type;
        computation_blocks_new_polymer[dep_v][dep_u].n_segment_compute = blocks[b].n_segment;
        computation_blocks_new_polymer[dep_v][dep_u].n_segment_offset = blocks[b].n_segment;
        computation_blocks_new_polymer[dep_v][dep_u].v_u.push_back(std::make_tuple(v,u));
        computation_blocks_new_polymer[dep_v][dep_u].n_repeated = computation_blocks_new_polymer[dep_v][dep_u].v_u.size();

        v_u_to_right_key[std::make_tuple(v,u)] = dep_u;
    }

    if (this->aggregate_propagator_computation)
    {

        // Aggregated keys
        std::map<std::string, std::vector<std::string>> aggregated_blocks;

        // Find aggregated branches in computation_blocks_new_polymer
        for(auto& item : computation_blocks_new_polymer)
        {
            // Left key and right keys
            std::string left_key = item.first;
            std::map<std::string, ComputationBlock> right_keys = item.second;

            // std::cout << "left_key: " << left_key << std::endl;

            // Aggregate propagators for given left key
            std::map<std::string, ComputationBlock> set_final;
            if (model_name == "continuous")
                set_final = PropagatorAnalyzer::aggregate_propagator_continuous_chain(right_keys);
            else if (model_name == "discrete")
                set_final = PropagatorAnalyzer::aggregate_propagator_discrete_chain(right_keys);
            else if (model_name == "")
                std::cout << "Chain model name is not set!" << std::endl;
            else
                std::cout << "Invalid model name: " << model_name << "!" << std::endl;

            // Replace the second map of computation_blocks_new_polymer with 'set_final'
            computation_blocks_new_polymer[left_key] = set_final;

            for(auto& item : set_final)
            {
                if(item.first[0] == '[' && item.second.n_segment_compute == item.second.n_segment_offset)
                    aggregated_blocks[left_key].push_back(item.first);
            }

            // Remove right keys from other left keys related to aggregated keys, and create new right keys

            substitute_right_keys(
                pc, v_u_to_right_key,
                computation_blocks_new_polymer,
                aggregated_blocks, left_key);
        }
    }

    // Add results to computation_blocks and computation_propagator_codes
    for(const auto& v_item : computation_blocks_new_polymer)
    {
        for(const auto& u_item : v_item.second)
        {
            std::string key_v = v_item.first;
            std::string key_u = u_item.first;

            int n_segment_compute = u_item.second.n_segment_compute;
            int n_segment_offset = u_item.second.n_segment_offset;
            int n_repeated = u_item.second.n_repeated;

            // Add blocks
            auto key = std::make_tuple(polymer_id, key_v, key_u);

            computation_blocks[key].monomer_type = PropagatorCode::get_monomer_type_from_key(key_v);
            computation_blocks[key].n_segment_compute = n_segment_compute;
            computation_blocks[key].n_segment_offset = n_segment_offset;
            computation_blocks[key].v_u = u_item.second.v_u;
            computation_blocks[key].n_repeated = n_repeated;

            // std::cout << "computation_blocks[key].n_repeated: " << key_v << ", " << key_u  << ", " << computation_blocks[key].n_repeated << std::endl;

            // Update propagators
            update_computation_propagator_map(computation_propagator_codes, key_v, n_segment_offset);
            update_computation_propagator_map(computation_propagator_codes, key_u, n_segment_compute);
        }
    }
}
std::map<std::string, ComputationBlock> PropagatorAnalyzer::aggregate_propagator_continuous_chain(std::map<std::string, ComputationBlock> not_aggregated_right_keys)
{
    // Example)
    // 0, B:
    //   N_offset, N_compute, R, C_R, 
    //   6, 6, 1, (C)B, 
    //   4, 4, 3, (D)B, 
    //   4, 4, 2, (E)B, 
    //   2, 2, 1, (F)B, 
    //
    //      ↓   Aggregation
    //  
    //   6, 6, 1, (C)B, 
    //   4, 0, 3, (D)B, 
    //   4, 0, 2, (E)B, 
    //   4, 4, 1, [(D)B0:3,(E)B0:2]B, 
    //   2, 2, 1, (F)B, 

    return aggregate_propagator_common(not_aggregated_right_keys, 0);

}
std::map<std::string, ComputationBlock> PropagatorAnalyzer::aggregate_propagator_discrete_chain(std::map<std::string, ComputationBlock> not_aggregated_right_keys)
{
    // Example)
    // 0, B:
    //   N_offset, N_compute, R, C_R, 
    //   6, 6, 1, (C)B, 
    //   4, 4, 3, (D)B, 
    //   4, 4, 2, (E)B, 
    //   2, 2, 1, (F)B, 
    //
    //      ↓   Aggregation
    //  
    //   6, 6, 1, (C)B, 
    //   4, 1, 3, (D)B, 
    //   4, 1, 2, (E)B, 
    //   4, 3, 1, [(D)B1:3,(E)B1:2]B, 
    //   2, 2, 1, (F)B, 

    return aggregate_propagator_common(not_aggregated_right_keys, 1);
}

std::map<std::string, ComputationBlock> PropagatorAnalyzer::aggregate_propagator_common(std::map<std::string, ComputationBlock> set_I, int minimum_n_segment)
{
    #ifndef NDEBUG
    std::cout << "--------- PropagatorAnalyzer::aggregate_propagator (before) -----------" << std::endl;
    std::cout << "--------- map ------------" << std::endl;
    for(const auto& item : set_I)
    {
        std::cout << item.second.n_segment_compute << ", " <<
                     item.first << ", " <<
                     item.second.n_segment_offset << ", ";
        for(const auto& v_u : item.second.v_u)
        {
            std::cout << "("
            + std::to_string(std::get<1>(v_u)) + ","
            + std::to_string(std::get<0>(v_u)) + ")" + ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "-----------------------" << std::endl;
    #endif

    std::set<int> set_n_compute;
    for(const auto& item: set_I)
        set_n_compute.insert(item.second.n_segment_compute);

    for(const int n_segment_current: set_n_compute)
    {
        // Add elements into set_S
        std::map<std::string, ComputationBlock, ComparePropagatorKey> set_S;
        for(const auto& item: set_I)
        {
            if (item.second.n_segment_compute == n_segment_current)
                set_S[item.first] = item.second;
        }

        // Skip if nothing to aggregate
        if (set_S.size() == 1)
            continue;

        // Update 'n_segment_compute'
        for(const auto& item: set_S)
            set_I[item.first].n_segment_compute = minimum_n_segment;
        
        // New 'n_segment_compute' and 'n_segment_offset'
        int n_segment_compute = n_segment_current - minimum_n_segment;
        int n_segment_offset  = n_segment_current - minimum_n_segment;
           
        // New 'v_u' and propagator key
        std::vector<std::tuple<int ,int>> v_u;
        std::string propagator_code = "[";
        bool is_first_sub_propagator = true;
        std::string dep_key;
        for(auto it = set_S.rbegin(); it != set_S.rend(); it++)
        {
            dep_key = it->first;

            // Update propagator key
            if(!is_first_sub_propagator)
                propagator_code += ",";
            else
                is_first_sub_propagator = false;
            propagator_code += dep_key + std::to_string(set_I[dep_key].n_segment_compute);

            // The number of repeats
            if (set_I[dep_key].n_repeated > 1)
                propagator_code += ":" + std::to_string(set_I[dep_key].n_repeated);

            // Compute the union of v_u
            std::vector<std::tuple<int ,int>> dep_v_u = set_I[dep_key].v_u;
            v_u.insert(v_u.end(), dep_v_u.begin(), dep_v_u.end());
        }
        std::string monomer_type = PropagatorCode::get_monomer_type_from_key(dep_key);
        propagator_code += "]" + monomer_type;

        // Add new aggregated key to set_I
        set_I[propagator_code].monomer_type = monomer_type;
        set_I[propagator_code].n_segment_compute = n_segment_compute;
        set_I[propagator_code].n_segment_offset = n_segment_offset;
        set_I[propagator_code].v_u = v_u;
        set_I[propagator_code].n_repeated = 1;

        #ifndef NDEBUG
        std::cout << "---------- PropagatorAnalyzer::aggregate_propagator (in progress) -----------" << std::endl;
        std::cout << "--------- map (" + std::to_string(set_I.size()) + ") -----------" << std::endl;
        for(const auto& item : set_I)
        {
            std::cout << item.second.n_segment_compute << ", " <<
                        item.first << ", " <<
                        item.second.n_segment_offset << ", ";
            for(const auto& v_u : item.second.v_u)
            {
                std::cout << "("
                + std::to_string(std::get<1>(v_u)) + ","
                + std::to_string(std::get<0>(v_u)) + ")" + ", ";
            }
            std::cout << std::endl;
        }
        std::cout << "-----------------------" << std::endl;
        #endif
    }
    return set_I;
}
bool PropagatorAnalyzer::is_aggregated() const
{
    return aggregate_propagator_computation;
}

void PropagatorAnalyzer::substitute_right_keys(
    Polymer& pc, 
    std::map<std::tuple<int, int>, std::string>& v_u_to_right_key,
    std::map<std::string, std::map<std::string, ComputationBlock>> & computation_blocks_new_polymer,
    std::map<std::string, std::vector<std::string>>& aggregated_blocks,
    std::string left_key)
{

    for(auto& aggregated_key : aggregated_blocks[left_key])
    {
        auto& computation_block = computation_blocks_new_polymer[left_key][aggregated_key];

        // std::cout << "aggregated_key: " << aggregated_key << std::endl;
        // For each v_u
        for(auto& old_v_u : computation_block.v_u)
        {
            // (old_u) -----  (old_v, u) ----- (v) 
            int old_v = std::get<0>(old_v_u);
            int old_u = std::get<1>(old_v_u);

            auto& vec_new_v = pc.get_adjacent_nodes()[old_v];
            // Remove 'old_u' from 'vec_new_v'
            vec_new_v.erase(std::remove(vec_new_v.begin(), vec_new_v.end(), old_u), vec_new_v.end());
            // std::cout << "(old_v_u): " << v <<  ", " << old_v << std::endl;

            // For each v_adj_node
            for(auto& v : vec_new_v)
            {
                int u = old_v;
                std::string dep_v = pc.get_propagator_key(v, u);
                std::string dep_u = pc.get_propagator_key(u, v);
                // std::cout << dep_v << ", " << dep_u << ", " << pc.get_block(v,u).n_segment << std::endl;

                // Make new key
                std::string new_u_key = "(" + aggregated_key
                    + std::to_string(computation_block.n_segment_compute);
                std::vector<std::string> sub_keys;

                for(auto& v_adj_node_dep : vec_new_v)
                {
                    if (v_adj_node_dep != v)
                        sub_keys.push_back(pc.get_propagator_key(v_adj_node_dep,u) + std::to_string(pc.get_block(v_adj_node_dep,u).n_segment));
                }
                std::sort(sub_keys.begin(),sub_keys.end());
                for(auto& item : sub_keys)
                    new_u_key += item;
                new_u_key += ")" + pc.get_block(v,u).monomer_type;

                // Remove 'v_u' from 'computation_blocks_new_polymer'
                // std::cout << "v_u_to_right_key[std::make_tuple(v,u)]: " << v_u_to_right_key[std::make_tuple(v,u)] << std::endl; 
                computation_blocks_new_polymer[dep_v].erase(v_u_to_right_key[std::make_tuple(v,u)]);

                // Add new key
                if (computation_blocks_new_polymer[dep_v].find(new_u_key) == computation_blocks_new_polymer[dep_v].end())
                {
                    computation_blocks_new_polymer[dep_v][new_u_key].monomer_type = pc.get_block(v,u).monomer_type;
                    computation_blocks_new_polymer[dep_v][new_u_key].n_segment_compute = pc.get_block(v,u).n_segment;
                    computation_blocks_new_polymer[dep_v][new_u_key].n_segment_offset = pc.get_block(v,u).n_segment;
                    computation_blocks_new_polymer[dep_v][new_u_key].v_u.push_back(std::make_tuple(v,u));

                    if (aggregated_key[0] == '[')
                        computation_blocks_new_polymer[dep_v][new_u_key].n_repeated = 1;
                    else
                        computation_blocks_new_polymer[dep_v][new_u_key].n_repeated = computation_block.n_repeated;

                    aggregated_blocks[dep_v].push_back(new_u_key);
                }
                else
                {
                    computation_blocks_new_polymer[dep_v][new_u_key].v_u.push_back(std::make_tuple(v,u));
                    int u0 = std::get<1>(computation_blocks_new_polymer[dep_v][new_u_key].v_u[0]);
                    if(u0 == u)
                        computation_blocks_new_polymer[dep_v][new_u_key].n_repeated += computation_block.n_repeated;
                }
                // std::cout << "dep_v, new_u_key, n_segment_compute, n_segment_offset : " << dep_v << ", " << new_u_key << ", " << n_segment_compute << ", " << n_segment_offset << std::endl;
            }
        }
    }
}

void PropagatorAnalyzer::update_computation_propagator_map(std::map<std::string, ComputationEdge, ComparePropagatorKey>& computation_propagator_codes, std::string new_key, int new_n_segment)
{
    if (computation_propagator_codes.find(new_key) == computation_propagator_codes.end())
    {
        computation_propagator_codes[new_key].deps = PropagatorCode::get_deps_from_key(new_key);
        computation_propagator_codes[new_key].monomer_type = PropagatorCode::get_monomer_type_from_key(new_key);
        computation_propagator_codes[new_key].max_n_segment = new_n_segment;
        computation_propagator_codes[new_key].height = PropagatorCode::get_height_from_key(new_key);
    }
    else
    {
        if (computation_propagator_codes[new_key].max_n_segment < new_n_segment)
            computation_propagator_codes[new_key].max_n_segment = new_n_segment;
    }
}
int PropagatorAnalyzer::get_n_computation_propagator_codes() const
{
    return computation_propagator_codes.size();
}
std::map<std::string, ComputationEdge, ComparePropagatorKey>& PropagatorAnalyzer::get_computation_propagator_codes()
{
    return computation_propagator_codes;
}
ComputationEdge& PropagatorAnalyzer::get_computation_propagator_code(std::string key)
{
    if (computation_propagator_codes.find(key) == computation_propagator_codes.end())
        throw_with_line_number("There is no such key (" + key + ").");

    return computation_propagator_codes[key];
}
std::map<std::tuple<int, std::string, std::string>, ComputationBlock>& PropagatorAnalyzer::get_computation_blocks()
{
    return computation_blocks;
}
ComputationBlock& PropagatorAnalyzer::get_computation_block(std::tuple<int, std::string, std::string> key)
{
    if (computation_blocks.find(key) == computation_blocks.end())
        throw_with_line_number("There is no such key (" + std::to_string(std::get<0>(key)) + ", " + 
            std::get<1>(key) + ", " + std::get<2>(key) + ").");

    return computation_blocks[key];
}
void PropagatorAnalyzer::display_blocks() const
{
    // Print blocks
    std::cout << "--------- Blocks ---------" << std::endl;
    std::cout << "Polymer id, left key:\n\taggregated, n_segment (offset, compute), right key, n_repeat, {v, u} list" << std::endl;

    const int MAX_PRINT_LENGTH = 500;
    std::tuple<int, std::string> v_tuple = std::make_tuple(-1, "");

    for(const auto& item : computation_blocks)
    {
        // Print polymer id, key1
        const std::string v_string = std::get<1>(item.first);
        if (v_tuple != std::make_tuple(std::get<0>(item.first), v_string))
        {
            std::cout << std::endl << std::to_string(std::get<0>(item.first)) + ", ";
            if (v_string.size() <= MAX_PRINT_LENGTH)
                std::cout << v_string;
            else
                std::cout << v_string.substr(0,MAX_PRINT_LENGTH-5) + " ... <omitted>, " ;
            std::cout << ":" << std::endl;
            v_tuple = std::make_tuple(std::get<0>(item.first), v_string);
        }

        // Print if aggregated
        const std::string u_string = std::get<2>(item.first);
        std::cout << "\t ";
        if (u_string.find('[') == std::string::npos)
            std::cout << "X, ";
        else
            std::cout << "O, ";
        // Print n_segment (offset, compute)
        std::cout << "(" + std::to_string(item.second.n_segment_offset) + ", " + std::to_string(item.second.n_segment_compute) + "), ";

        // Print key2
        if (u_string.size() <= MAX_PRINT_LENGTH)
            std::cout << u_string;
        else
            std::cout << u_string.substr(0,MAX_PRINT_LENGTH-5) + " ... <omitted>" ;

        // Print n_repeat
        std::cout << ", " + std::to_string(item.second.n_repeated);        

        // Print v_u list
        for(const auto& v_u : item.second.v_u)
        {
            std::cout << ", {"
            + std::to_string(std::get<1>(v_u)) + ","
            + std::to_string(std::get<0>(v_u)) + "}";
        }
        std::cout << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;
}
void PropagatorAnalyzer::display_propagators() const
{
    // Print propagators
    std::vector<std::tuple<std::string, int, int>> sub_deps;
    int total_segments = 0;

    std::cout << "--------- Propagators ---------" << std::endl;
    std::cout << "Key:\n\taggregated, max_n_segment, height" << std::endl;
    
    for(const auto& item : computation_propagator_codes)
    {
        total_segments += item.second.max_n_segment;

        const int MAX_PRINT_LENGTH = 500;

        if (item.first.size() <= MAX_PRINT_LENGTH)
            std::cout << item.first;
        else
            std::cout << item.first.substr(0,MAX_PRINT_LENGTH-5) + " ... <omitted> " ;

        std::cout << ":\n\t ";
        if (item.first.find('[') == std::string::npos)
            std::cout << "X, ";
        else
            std::cout << "O, ";
        std::cout << item.second.max_n_segment << ", " << item.second.height << std::endl;
    }
    std::cout << "Total number of modified diffusion equation (or integral equation for discrete chain model) steps to compute propagators: " << total_segments << std::endl;
    std::cout << "------------------------------------" << std::endl;
}

void PropagatorAnalyzer::display_sub_propagators() const
{
    // Print sub propagators
    std::vector<std::tuple<std::string, int, int>> sub_deps;
    int total_segments = 0;
    std::cout << "--------- Propagators ---------" << std::endl;
    std::cout << "Key:\n\taggregated, max_n_segment, height, deps," << std::endl;
    
    for(const auto& item : computation_propagator_codes)
    {
        total_segments += item.second.max_n_segment;

        std::cout << item.first;
        std::cout << ":\n\t ";
        if (item.first.find('[') == std::string::npos)
            std::cout << "X, ";
        else
            std::cout << "O, ";
        std::cout << item.second.max_n_segment << ", " << item.second.height;

        sub_deps = PropagatorCode::get_deps_from_key(item.first);
        for(size_t i=0; i<sub_deps.size(); i++)
        {
            std::cout << ", "  << std::get<0>(sub_deps[i]) << ":" << std::get<1>(sub_deps[i]);
        }
        std::cout << std::endl;
    }
    std::cout << "Total number of modified diffusion equation (or integral equation for discrete chain model) steps to compute propagators: " << total_segments << std::endl;
    std::cout << "------------------------------------" << std::endl;
}