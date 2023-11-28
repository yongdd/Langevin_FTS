#include <iostream>
#include <cctype>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stack>
#include <set>

#include "Molecules.h"
#include "Exception.h"

//----------------- Constructor ----------------------------
Molecules::Molecules(
    std::string model_name, double ds, std::map<std::string, double> bond_lengths, bool reduce_propagator_computation)
{
    // Checking chain model
    std::transform(model_name.begin(), model_name.end(), model_name.begin(),
                   [](unsigned char c)
    {
        return std::tolower(c);
    });

    if (model_name != "continuous" && model_name != "discrete")
    {
        throw_with_line_number(model_name + " is an invalid chain model. This must be 'Continuous' or 'Discrete'.");
    }
    this->model_name = model_name;

    // Save variables
    try
    {
        this->ds = ds;
        this->bond_lengths = bond_lengths;
        this->reduce_propagator_computation = reduce_propagator_computation;
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }
}
void Molecules::add_polymer(
    double volume_fraction,
    std::vector<BlockInput> block_inputs,
    std::map<int, std::string> chain_end_to_q_init)
{
    std::string propagator_code;
    distinct_polymers.push_back(Polymer(ds, bond_lengths, 
        volume_fraction, block_inputs, chain_end_to_q_init));

    Polymer& pc = distinct_polymers.back();

    // Construct starting vertices 'v', ending vertices 'u', 
    std::vector<int> v;
    std::vector<int> u;
    for(size_t i=0; i<block_inputs.size(); i++)
    {
        v.push_back(block_inputs[i].v);
        u.push_back(block_inputs[i].u);
    }

    // Generate propagator code for each block and each direction
    std::map<std::pair<int, int>, std::pair<std::string, int>> memory;
    for (int i=0; i<pc.get_n_blocks(); i++)
    {
        propagator_code = generate_propagator_code(
            memory,
            pc.get_blocks(),
            pc.get_adjacent_nodes(),
            pc.get_block_indexes(),
            chain_end_to_q_init,
            v[i], u[i]).first;
        pc.set_propagator_key(propagator_code, v[i], u[i]);

        propagator_code = generate_propagator_code(
            memory,
            pc.get_blocks(),
            pc.get_adjacent_nodes(),
            pc.get_block_indexes(),
            chain_end_to_q_init,
            u[i], v[i]).first;
        pc.set_propagator_key(propagator_code, u[i], v[i]);
    }

    // Temporary map for the new polymer
    std::map<std::tuple<int, std::string>, std::map<std::string, EssentialBlock >> essential_blocks_new_polymer;

    // Find essential_blocks in new_polymer
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

        auto key1 = std::make_tuple(distinct_polymers.size()-1, dep_v);
        essential_blocks_new_polymer[key1][dep_u].monomer_type = blocks[b].monomer_type;
        essential_blocks_new_polymer[key1][dep_u].n_segment_allocated = blocks[b].n_segment;
        essential_blocks_new_polymer[key1][dep_u].n_segment_offset    = 0;
        essential_blocks_new_polymer[key1][dep_u].n_segment_original  = blocks[b].n_segment;
        essential_blocks_new_polymer[key1][dep_u].v_u.push_back(std::make_tuple(v,u));
    }

    if (this->reduce_propagator_computation)
    {
        // Find superposed branches in essential_blocks_new_polymer
        std::map<std::tuple<int, std::string>, std::map<std::string, EssentialBlock >> superposed_blocks;
        for(auto& item : essential_blocks_new_polymer)
        {
            
            std::vector<std::tuple<int, int>> total_v_u_list;
            // Find all (v,u) pairs in superposed_blocks for the given key
            for(auto& second_key : superposed_blocks[item.first]) // map <tuple, v_u_vector>
            {
                for(auto& superposition_v_u : second_key.second.v_u) 
                {
                    total_v_u_list.push_back(superposition_v_u);
                    // std::cout << "(v_u): " << std::get<0>(superposition_v_u) <<  ", " << std::get<1>(superposition_v_u) << std::endl;
                }
            }

            // Remove keys of second map in essential_blocks_new_polymer, which exist in superposed_blocks.
            for(auto it = item.second.cbegin(); it != item.second.cend();) // map <tuple, v_u_vector>
            {
                bool removed = false;
                for(auto& v_u : total_v_u_list)
                {
                    if ( std::find(it->second.v_u.begin(), it->second.v_u.end(), v_u) != it->second.v_u.end())
                    {
                        it = item.second.erase(it);
                        removed = true;
                        break;
                    }
                }
                if (!removed)
                    ++it;
            }

            // After the removal is done, add the superposed branches
            for(auto& second_key : superposed_blocks[item.first]) 
                essential_blocks_new_polymer[item.first][second_key.first] = second_key.second;

            // Superpose propagators for given key
            // If the number of elements in the second map is only 1, it will return the map without superposition.
            // Not all elements of superposed_second_map are superposed.
            std::map<std::string, EssentialBlock> superposed_second_map;
            if (model_name == "continuous")
                superposed_second_map = superpose_propagator_of_continuous_chain(item.second);
            else if (model_name == "discrete")
                superposed_second_map = superpose_propagator_of_discrete_chain(item.second);

            // Replace the second map of essential_blocks_new_polymer with superposed_second_map
            essential_blocks_new_polymer[item.first].clear();
            for(auto& superposed_propagator_code : superposed_second_map)
                essential_blocks_new_polymer[item.first][superposed_propagator_code.first] = superposed_propagator_code.second;

            // For each superposed_propagator_code
            for(auto& superposed_propagator_code : superposed_second_map)
            {
                int n_segment_allocated = superposed_propagator_code.second.n_segment_allocated;
                std::string dep_key = superposed_propagator_code.first;
                int n_segment_offset = superposed_propagator_code.second.n_segment_offset;
                int n_segment_original = superposed_propagator_code.second.n_segment_original;

                // Skip, if it is not superposed
                if ( dep_key[0] != '[' || n_segment_offset+n_segment_allocated != n_segment_original)
                    continue;

                // For each v_u 
                for(auto& v_u : superposed_propagator_code.second.v_u)
                {
                    auto& v_adj_nodes = pc.get_adjacent_nodes()[std::get<0>(v_u)];
                    // For each v_adj_node
                    for(auto& v_adj_node : v_adj_nodes)
                    {
                        if (v_adj_node != std::get<1>(v_u))
                        {
                            // std::cout << "(v_u): " << v_adj_node <<  ", " << std::get<0>(v_u) << std::endl;
                            int v = v_adj_node;
                            int u = std::get<0>(v_u);

                            std::string dep_v = pc.get_propagator_key(v, u);
                            std::string dep_u = pc.get_propagator_key(u, v);
                            // std::cout << dep_v << ", " << dep_u << ", " << pc.get_block(v,u).n_segment << std::endl;

                            auto key = std::make_tuple(distinct_polymers.size()-1, dep_v);
                            // pc.get_block(v,u).monomer_type
                            // pc.get_block(v,u).n_segment

                            // Make new key
                            std::string new_u_key = "(" + dep_key
                                + std::to_string(n_segment_allocated);

                            for(auto& v_adj_node_dep : v_adj_nodes)
                            {
                                if (v_adj_node_dep != v && v_adj_node_dep != std::get<1>(v_u))
                                    new_u_key += pc.get_block(v_adj_node_dep,u).monomer_type + std::to_string(pc.get_block(v_adj_node_dep,u).n_segment);
                                
                            }
                            new_u_key += ")" + pc.get_block(v,u).monomer_type;

                            // Add the new key
                            superposed_blocks[key][new_u_key].monomer_type = pc.get_block(v,u).monomer_type;
                            superposed_blocks[key][new_u_key].n_segment_allocated = pc.get_block(v,u).n_segment;
                            superposed_blocks[key][new_u_key].n_segment_offset = 0;
                            superposed_blocks[key][new_u_key].n_segment_original = pc.get_block(v,u).n_segment;
                            superposed_blocks[key][new_u_key].v_u.push_back(std::make_tuple(v,u));
                        }
                    }
                }
            }
        }
    }

    // Add results to essential_blocks and essential_propagator_codes
    for(const auto& v_item : essential_blocks_new_polymer)
    {
        for(const auto& u_item : v_item.second)
        {
            int polymer_id = std::get<0>(v_item.first);

            std::string key_v = std::get<1>(v_item.first);
            std::string key_u = u_item.first;

            int n_segment_allocated = u_item.second.n_segment_allocated;
            int n_segment_offset = u_item.second.n_segment_offset;
            int n_segment_original = u_item.second.n_segment_original;

            // Add blocks
            auto key = std::make_tuple(polymer_id, key_v, key_u);

            essential_blocks[key].monomer_type = Molecules::get_monomer_type_from_key(key_v);
            essential_blocks[key].n_segment_allocated = n_segment_allocated;
            essential_blocks[key].n_segment_offset    = n_segment_offset;
            essential_blocks[key].n_segment_original  = n_segment_original;
            essential_blocks[key].v_u = u_item.second.v_u;

            // Update propagators
            update_essential_propagator_code(essential_propagator_codes, key_v, n_segment_original);
            update_essential_propagator_code(essential_propagator_codes, key_u, n_segment_allocated);
        }
    }
}
std::string Molecules::get_model_name() const
{
    return model_name;
}
double Molecules::get_ds() const
{
    return ds;
}
bool Molecules::is_using_superposition() const
{
    return reduce_propagator_computation;
}
int Molecules::get_n_polymer_types() const
{
    return distinct_polymers.size();
}
Polymer& Molecules::get_polymer(const int p)
{
    return distinct_polymers[p];
}
const std::map<std::string, double>& Molecules::get_bond_lengths() const
{
    return bond_lengths;
}
std::pair<std::string, int> Molecules::generate_propagator_code(
    std::map<std::pair<int, int>, std::pair<std::string, int>>& memory,
    std::vector<Block>& blocks,
    std::map<int, std::vector<int>>& adjacent_nodes,
    std::map<std::pair<int, int>, int>& edge_to_block_index,
    std::map<int, std::string>& chain_end_to_q_init,
    int in_node, int out_node)
{
    std::vector<std::string> edge_text;
    std::vector<std::pair<std::string,int>> edge_dict;
    std::pair<std::string,int> text_and_segments;

    // If it is already computed
    if (memory.find(std::make_pair(in_node, out_node)) != memory.end())
        return memory[std::make_pair(in_node, out_node)];

    // Explore child blocks
    //std::cout << "[" + std::to_string(in_node) + ", " +  std::to_string(out_node) + "]:";
    for(size_t i=0; i<adjacent_nodes[in_node].size(); i++)
    {
        if (adjacent_nodes[in_node][i] != out_node)
        {
            //std::cout << "(" << in_node << ", " << adjacent_nodes[in_node][i] << ")";
            auto v_u_pair = std::make_pair(adjacent_nodes[in_node][i], in_node);
            if (memory.find(v_u_pair) != memory.end())
                text_and_segments = memory[v_u_pair];
            else
            {
                text_and_segments = generate_propagator_code(
                    memory, blocks, adjacent_nodes, edge_to_block_index,
                    chain_end_to_q_init,
                    adjacent_nodes[in_node][i], in_node);
                memory[v_u_pair] = text_and_segments;
            }
            edge_text.push_back(text_and_segments.first + std::to_string(text_and_segments.second));
            edge_dict.push_back(text_and_segments);
            //std::cout << text_and_segments.first << " " << text_and_segments.second << std::endl;
        }
    }

    // Merge code of child blocks
    std::string text;
    if(edge_text.size() == 0)
    {
        // If in_node does not exist in chain_end_to_q_init
        if (chain_end_to_q_init.find(in_node) == chain_end_to_q_init.end())
            text = "";

        else
        {
            text = "{" + chain_end_to_q_init[in_node] + "}";
        }
    }
    else
    {
        std::sort(edge_text.begin(), edge_text.end());
        text += "(";
        for(size_t i=0; i<edge_text.size(); i++)
            text += edge_text[i];
        text += ")";
    }

    // Add monomer_type at the end of text code
    text += blocks[edge_to_block_index[std::make_pair(in_node, out_node)]].monomer_type;
    auto text_and_segments_total = std::make_pair(text, blocks[edge_to_block_index[std::make_pair(in_node, out_node)]].n_segment);
    memory[std::make_pair(in_node, out_node)] = text_and_segments_total;

    return text_and_segments_total;
}
void Molecules::update_essential_propagator_code(std::map<std::string, EssentialEdge, ComparePropagatorKey>& essential_propagator_codes, std::string new_key, int new_n_segment)
{
    if (essential_propagator_codes.find(new_key) == essential_propagator_codes.end())
    {
        essential_propagator_codes[new_key].deps = Molecules::get_deps_from_key(new_key);
        essential_propagator_codes[new_key].monomer_type = Molecules::get_monomer_type_from_key(new_key);
        essential_propagator_codes[new_key].max_n_segment = new_n_segment;
        essential_propagator_codes[new_key].height = Molecules::get_height_from_key(new_key);
    }
    else
    {
        if (essential_propagator_codes[new_key].max_n_segment < new_n_segment)
            essential_propagator_codes[new_key].max_n_segment = new_n_segment;
    }
}
std::map<std::string, EssentialBlock> Molecules::superpose_propagator_of_continuous_chain(std::map<std::string, EssentialBlock> not_superposed_yet_second_map)
{
    // Example)
    // 0, B:
    //   6, 0, 6, (C)B, 1,
    //   4, 0, 4, (D)B, 3,
    //   4, 0, 4, (E)B, 2,
    //   2, 0, 2, (F)B, 1,
    //
    //      ↓   Superposition
    //  
    //   6, 0, 2, (C)B, 1,  // done
    //   4, 0, 0, (D)B, 3,  // done
    //   4, 0, 0, (E)B, 2,  // done
    //   2, 0, 2, (F)B, 1,
    //   6, 2, 4, [(C)B2:1,(D)B0:3,(E)B0:2]B,
    //
    //      ↓   Superposition
    //  
    //   6, 0, 2, (C)B, 1,  // done
    //   4, 0, 0, (D)B, 3,  // done
    //   4, 0, 0, (E)B, 2,  // done
    //   2, 0, 0, (F)B, 1,  // done
    //   6, 2, 2, [(C)B2:1,(D)B0:3,(E)B0:2]B,             // done
    //   6, 4, 2, [[(C)B2:1,(D)B0:3,(E)B0:2]B2,(F)B2:1]B  // done

    std::map<std::string, EssentialBlock> remaining_keys;
    std::map<std::string, EssentialBlock> superposed_second_map;
    std::map<std::string, EssentialBlock> superposed_second_map_total;

    // Because of our SimpsonRule implementation, whose weights of odd number n_segments and even number n_segments are slightly different,
    // superpositions for blocks of odd number and of even number are separately performed.

    // For even number
    for(const auto& item : not_superposed_yet_second_map)
    {
        if (item.second.n_segment_allocated % 2 == 0)
            remaining_keys[item.first] = item.second;
    }
    superposed_second_map_total = superpose_propagator_common(remaining_keys, 0);

    // For odd number
    remaining_keys.clear();
    for(const auto& item : not_superposed_yet_second_map)
    {
        if (item.second.n_segment_allocated % 2 == 1)
            remaining_keys[item.first] = item.second;
    }
    superposed_second_map = superpose_propagator_common(remaining_keys, 0);

    // Merge maps
    superposed_second_map_total.insert(std::begin(superposed_second_map), std::end(superposed_second_map));

    return superposed_second_map_total;

}
std::map<std::string, EssentialBlock> Molecules::superpose_propagator_of_discrete_chain(std::map<std::string, EssentialBlock> not_superposed_yet_second_map)
{

    // Example)
    // 0, B:
    //   6, 0, 6, (C)B, 1,
    //   4, 0, 4, (D)B, 3,
    //   4, 0, 4, (E)B, 2,
    //   2, 0, 2, (F)B, 1,
    //
    //      ↓   Superposition
    //  
    //   6, 0, 3, (C)B, 1,  // done
    //   4, 0, 1, (D)B, 3,  // done
    //   4, 0, 1, (E)B, 2,  // done
    //   2, 0, 2, (F)B, 1,
    //   6, 3, 3, [(C)B3:1,(D)B1:3,(E)B1:2]B,
    //
    //      ↓   Superposition
    //  
    //   6, 0, 3, (C)B, 1,  // done
    //   4, 0, 1, (D)B, 3,  // done
    //   4, 0, 1, (E)B, 2,  // done
    //   2, 0, 1, (F)B, 1,  // done
    //   6, 3, 2, [(C)B3:1,(D)B1:3,(E)B1:2]B,             // done
    //   6, 5, 1, [[(C)B3:1,(D)B1:3,(E)B1:2]B2,(F)B1:1]B  // done

    // std::map<std::string, EssentialBlock> remaining_keys;
    // for(const auto& item : not_superposed_yet_second_map)
    // {
    //     remaining_keys[item.first] = item.second;
    // }

    return superpose_propagator_common(not_superposed_yet_second_map, 1);
}

std::map<std::string, EssentialBlock> Molecules::superpose_propagator_common(std::map<std::string, EssentialBlock> remaining_keys, int minimum_n_segment)
{
    int current_n_segment;
    int n_segment_allocated;
    int n_segment_offset;
    int n_segment_original;

    std::string superposed_propagator_code;
    std::vector<std::tuple<int ,int>> v_u_total;

    std::map<std::string, EssentialBlock> superposed_second_map;
    // Tuple <n_segment_allocated, key, n_segment_offset, n_segment_original, v_u_list>
    std::vector<std::tuple<int, std::string, int, int, std::vector<std::tuple<int ,int>>>> same_superposition_level_list;

    // std::cout << "---------map------------" << std::endl;
    // for(const auto& item : remaining_keys)
    // {
    //     std::cout << item.second.n_segment_allocated << ", " <<
    //                  item.first << ", " <<
    //                  item.second.n_segment_offset << ", " <<
    //                  item.second.n_segment_original << ", ";
    //     for(const auto& v_u : item.second.v_u)
    //     {
    //         std::cout << "("
    //         + std::to_string(std::get<0>(v_u)) + ","
    //         + std::to_string(std::get<1>(v_u)) + ")" + ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "-----------------------" << std::endl;

    // Int count = 0;
    current_n_segment = 0;
    for(const auto& item: remaining_keys)
        current_n_segment = std::max(current_n_segment, item.second.n_segment_allocated);
    while(!remaining_keys.empty())
    {
        // Count ++;
        // If (count == 10)
        //     break;
        // std::cout << "------remaining_keys------" << std::endl;
        // for(const auto& item: remaining_keys)
        // {
        //     std::cout << item.second.n_segment_allocated << ", "  
        //               << item.first << ", " 
        //               << item.second.n_segment_offset << ", "
        //               << item.second.n_segment_original << ", " << std::endl;
        // }
        // std::cout << "-------------" << std::endl;

        // Find keys that have the same superposition level from remaining_keys
        std::set<int, std::greater<int>> n_segment_set; // for finding the largest n_segment that is not in same_superposition_level_list.
        for (auto it = remaining_keys.cbegin(); it != remaining_keys.cend(); )
        {
            bool erased = false;
            if (it->second.n_segment_allocated <= 1)
            {
                superposed_second_map[it->first] = it->second;
                remaining_keys.erase(it++);
                erased = true;
            }

            if (!erased)
            {
                if (current_n_segment <= it->second.n_segment_allocated)
                    same_superposition_level_list.push_back(std::make_tuple(
                        it->second.n_segment_allocated,
                        it->first,
                        it->second.n_segment_offset,
                        it->second.n_segment_original,
                        it->second.v_u));
                else
                    n_segment_set.insert(it->second.n_segment_allocated);
                ++it;
            }
        }
        // If it empty, decrease current_n_segment.
        if (same_superposition_level_list.empty())
            current_n_segment = *std::next(n_segment_set.begin(), 0);
        else
        {
            // std::cout << "------same_superposition_level_list------" << std::endl;
            // for(size_t i=0; i<same_superposition_level_list.size(); i++)
            //     std::cout << std::get<0>(same_superposition_level_list[i]) << ", " << std::get<1>(same_superposition_level_list[i]) << std::endl;
            // std::cout << "------------------------------------" << std::endl;

            v_u_total.clear();
            // If there is only one element
            if (same_superposition_level_list.size() == 1)
            {
                //  no the second largest element
                if (n_segment_set.size() == 0)
                {
                    // Add to map
                    superposed_second_map[std::get<1>(same_superposition_level_list[0])].monomer_type = Molecules::get_monomer_type_from_key(std::get<1>(same_superposition_level_list[0]));
                    superposed_second_map[std::get<1>(same_superposition_level_list[0])].n_segment_allocated = std::get<0>(same_superposition_level_list[0]);
                    superposed_second_map[std::get<1>(same_superposition_level_list[0])].n_segment_offset    = std::get<2>(same_superposition_level_list[0]);
                    superposed_second_map[std::get<1>(same_superposition_level_list[0])].n_segment_original  = std::get<3>(same_superposition_level_list[0]);
                    superposed_second_map[std::get<1>(same_superposition_level_list[0])].v_u                 = std::get<4>(same_superposition_level_list[0]);

                    // Erase element
                    remaining_keys.erase(std::get<1>(same_superposition_level_list[0]));
                }
                // Lower 'current_n_segment' to the next level and repeat
                else
                    current_n_segment = *std::next(n_segment_set.begin(), 0);
            }
            // Superposition
            else
            {
                // sort same_superposition_level_list with height in descending order
                std::sort(same_superposition_level_list.begin(), same_superposition_level_list.end(),
                    [](auto const &t1, auto const &t2)
                        {
                            return Molecules::get_height_from_key(std::get<1>(t1)) > Molecules::get_height_from_key(std::get<1>(t2));
                        }
                );

                // Add one by one
                std::string dep_key;
                std::vector<std::tuple<int ,int>> dep_v_u ;
                int n_segment_offset_max = 0;
                int n_segment_original_max = 0;

                for(size_t i=0; i<same_superposition_level_list.size(); i++)
                {
                    n_segment_allocated = std::get<0>(same_superposition_level_list[i]) - current_n_segment + minimum_n_segment;

                    dep_key = std::get<1>(same_superposition_level_list[i]);
                    n_segment_offset = std::get<2>(same_superposition_level_list[i]);
                    n_segment_original = std::get<3>(same_superposition_level_list[i]);
                    dep_v_u = std::get<4>(same_superposition_level_list[i]);

                    n_segment_offset_max = std::max(n_segment_offset_max, n_segment_offset + n_segment_allocated);
                    n_segment_original_max = std::max(n_segment_original_max, n_segment_original);

                    v_u_total.insert(v_u_total.end(), dep_v_u.begin(), dep_v_u.end());
                    if (i==0)
                        superposed_propagator_code = "[" + dep_key + std::to_string(n_segment_allocated);
                    else
                        superposed_propagator_code += "," + dep_key + std::to_string(n_segment_allocated);

                    if (dep_key.find('[') == std::string::npos)
                        superposed_propagator_code += ":" + std::to_string(dep_v_u.size());

                    // Add to map
                    superposed_second_map[std::get<1>(same_superposition_level_list[i])].monomer_type = Molecules::get_monomer_type_from_key(dep_key);
                    superposed_second_map[std::get<1>(same_superposition_level_list[i])].n_segment_allocated = n_segment_allocated;
                    superposed_second_map[std::get<1>(same_superposition_level_list[i])].n_segment_offset    = n_segment_offset;
                    superposed_second_map[std::get<1>(same_superposition_level_list[i])].n_segment_original  = n_segment_original;
                    superposed_second_map[std::get<1>(same_superposition_level_list[i])].v_u                 = dep_v_u;
                }
                superposed_propagator_code += "]" + Molecules::get_monomer_type_from_key(dep_key);
                n_segment_allocated = current_n_segment - minimum_n_segment;

                // Add to remaining_keys
                remaining_keys[superposed_propagator_code].monomer_type = Molecules::get_monomer_type_from_key(superposed_propagator_code);
                remaining_keys[superposed_propagator_code].n_segment_allocated = n_segment_allocated;
                remaining_keys[superposed_propagator_code].n_segment_offset    = n_segment_offset_max;
                remaining_keys[superposed_propagator_code].n_segment_original  = n_segment_original_max;
                remaining_keys[superposed_propagator_code].v_u                 = v_u_total;

                // Erase elements
                for(size_t i=0; i<same_superposition_level_list.size(); i++)
                    remaining_keys.erase(std::get<1>(same_superposition_level_list[i]));
            }
            same_superposition_level_list.clear();
        }
    }
    return superposed_second_map;
}

int Molecules::get_n_essential_propagator_codes() const
{
    return essential_propagator_codes.size();
}
std::vector<std::tuple<std::string, int, int>> Molecules::get_deps_from_key(std::string key)
{
    std::vector<std::tuple<std::string, int, int>> sub_deps;
    int sub_n_segment;
    std::string sub_key;
    int sub_n_repeated;

    bool is_reading_key = true;
    bool is_reading_n_segment = false;
    bool is_reading_n_repeated = false;

    int key_start = 1;
    int brace_count = 0;

    for(size_t i=0; i<key.size();i++)
    {
        // It was reading key and have found a digit
        if( isdigit(key[i]) && is_reading_key && brace_count == 1 )
        {
            // std::cout << "key_to_deps1" << std::endl;
            sub_key = key.substr(key_start, i-key_start);
            // std::cout << sub_key << "= " << key_start << ", " << i  << std::endl;

            is_reading_key = false;
            is_reading_n_segment = true;

            key_start = i;
        }
        // It was reading n_segment and have found a ':'
        else if( key[i]==':' && is_reading_n_segment && brace_count == 1 )
        {
            // std::cout << "key_to_deps2" << std::endl;
            sub_n_segment = std::stoi(key.substr(key_start, i-key_start));
            // std::cout << sub_key << "= " << key_start << ", " << i  << ", " << key.substr(key_start, i-key_start) << std::endl;

            is_reading_n_segment = false;
            is_reading_n_repeated = true;

            key_start = i+1;
        }
        // It was reading n_segment and have found a comma
        else if( key[i]==',' && is_reading_n_segment && brace_count == 1 )
        {
            // std::cout << "key_to_deps3" << std::endl;
            sub_n_segment = std::stoi(key.substr(key_start, i-key_start));
            // std::cout << sub_key << "= " << key_start << ", " << i  << ", " << key.substr(key_start, i-key_start) << std::endl;
            sub_deps.push_back(std::make_tuple(sub_key, sub_n_segment, 1));

            is_reading_n_segment = false;
            is_reading_key = true;

            key_start = i+1;
        }
        // It was reading n_repeated and have found a comma
        else if( key[i]==',' && is_reading_n_repeated && brace_count == 1 )
        {
            // std::cout << "key_to_deps4" << std::endl;
            sub_n_repeated = std::stoi(key.substr(key_start, i-key_start));
            // std::cout << sub_key << "= " << key_start << ", " << i  << ", " << key.substr(key_start, i-key_start) << std::endl;
            sub_deps.push_back(std::make_tuple(sub_key, sub_n_segment, sub_n_repeated));

            is_reading_n_repeated = false;
            is_reading_key = true;

            key_start = i+1;
        }
        // It was reading n_repeated and have found a non-digit
        else if( !isdigit(key[i]) && is_reading_n_repeated && brace_count == 1)
        {
            // std::cout << "key_to_deps5" << std::endl;
            sub_n_repeated = std::stoi(key.substr(key_start, i-key_start));
            // std::cout << sub_key << "= " << key_start << ", " << i  << ", " << key.substr(key_start, i-key_start) << std::endl;
            sub_deps.push_back(std::make_tuple(sub_key, sub_n_segment, sub_n_repeated));

            is_reading_n_repeated = false;
            is_reading_key = true;

            key_start = i;
        }
        // It was reading n_segment and have found a non-digit
        else if( !isdigit(key[i]) && is_reading_n_segment && brace_count == 1)
        {
            // std::cout << "key_to_deps6" << std::endl;
            sub_n_segment = std::stoi(key.substr(key_start, i-key_start));
            // std::cout << sub_key << "= " << key_start << ", " << i  << ", " << key.substr(key_start, i-key_start) << std::endl;
            sub_deps.push_back(std::make_tuple(sub_key, sub_n_segment, 1));

            is_reading_n_segment = false;
            is_reading_key = true;

            key_start = i;
        }
        if(key[i] == '(' || key[i] == '[')
            brace_count++;
        else if(key[i] == ')' || key[i] == ']')
            brace_count--;
    }
    return sub_deps;
}

std::string Molecules::remove_monomer_type_from_key(std::string key)
{
    if (key[0] != '[' && key[0] != '(' && key[0] != '{')
    {
        return "";
    }
    else
    {
        int brace_count = 0;
        int species_idx = 0;
        for(size_t i=0; i<key.size();i++)
        {
            if (key[i] == '[' || key[i] == '(' || key[i] == '{')
            {
                brace_count++;
            }
            else if (key[i] == ']' || key[i] == ')' || key[i] == '}')
            {
                brace_count--;
                if (brace_count == 0)
                {
                    species_idx=i;
                    break;
                }
            }
        }
        // std::cout << "key.substr(1, species_idx): " << key.substr(1, species_idx-1) << std::endl;
        return key.substr(1, species_idx-1);
    }
}

std::string Molecules::get_monomer_type_from_key(std::string key)
{
    int key_start = 0;
    for(int i=key.size()-1; i>=0;i--)
    {
        //std::cout << key[i] << std::endl;
        if(key[i] == ')' || key[i] == ']' || key[i] == '}')
        {
            key_start=i+1;
            break;
        }
    }
    //std::cout << key.substr(key_start, key.size()-key_start) << std::endl;
    return key.substr(key_start, key.size()-key_start);
}
std::string Molecules::get_q_input_idx_from_key(std::string key)
{
    if (key[0] != '{')
        throw_with_line_number("There is no related initial condition in key (" + key + ").");

    int key_start = 0;
    for(int i=key.size()-1; i>=0;i--)
    {
        if(key[i] == '}')
        {
            key_start=i;
            break;
        }
    }
    // std::cout << key.substr(1, key_start-1) << std::endl;
    return key.substr(1, key_start-1);
}
int Molecules::get_height_from_key(std::string key)
{
    int height_count = 0;
    for(size_t i=0; i<key.size();i++)
    {
        if (key[i] == '[' || key[i] == '(')
            height_count++;
        else
            break;
    }
    return height_count;
}
std::map<std::string, EssentialEdge, ComparePropagatorKey>& Molecules::get_essential_propagator_codes()
{
    return essential_propagator_codes;
}
EssentialEdge& Molecules::get_essential_propagator_code(std::string key)
{
    if (essential_propagator_codes.find(key) == essential_propagator_codes.end())
        throw_with_line_number("There is no such key (" + key + ").");

    return essential_propagator_codes[key];
}
std::map<std::tuple<int, std::string, std::string>, EssentialBlock>& Molecules::get_essential_blocks()
{
    return essential_blocks;
}
EssentialBlock& Molecules::get_essential_block(std::tuple<int, std::string, std::string> key)
{
    if (essential_blocks.find(key) == essential_blocks.end())
        throw_with_line_number("There is no such key (" + std::to_string(std::get<0>(key)) + ", " + 
            std::get<1>(key) + ", " + std::get<2>(key) + ").");

    return essential_blocks[key];
}
void Molecules::display_blocks() const
{
    // Print blocks
    std::cout << "--------- Blocks ---------" << std::endl;
    std::cout << "Polymer id, key1:\n\tsuperposed, n_segment (original, offset, allocated), key2, {v, u} list" << std::endl;

    const int MAX_PRINT_LENGTH = 500;
    std::tuple<int, std::string> v_tuple = std::make_tuple(-1, "");

    for(const auto& item : essential_blocks)
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

        // Print if superposed
        const std::string u_string = std::get<2>(item.first);
        std::cout << "\t ";
        if (u_string.find('[') == std::string::npos)
            std::cout << "X, ";
        else
            std::cout << "O, ";
        // Print n_segment (original, offset, allocated)
        std::cout << "(" + std::to_string(item.second.n_segment_original) + ", "+ std::to_string(item.second.n_segment_offset) + ", " + std::to_string(item.second.n_segment_allocated) + "), ";

        // Print key2
        if (u_string.size() <= MAX_PRINT_LENGTH)
            std::cout << u_string;
        else
            std::cout << u_string.substr(0,MAX_PRINT_LENGTH-5) + " ... <omitted>" ;

        // Print v_u list
        for(const auto& v_u : item.second.v_u)
        {
            std::cout << ", {"
            + std::to_string(std::get<0>(v_u)) + ","
            + std::to_string(std::get<1>(v_u)) + "}";
        }
        std::cout << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;
}
void Molecules::display_propagators() const
{
    // Print propagators
    std::vector<std::tuple<std::string, int, int>> sub_deps;
    int total_segments = 0;

    std::cout << "--------- Propagators ---------" << std::endl;
    std::cout << "Key:\n\tsuperposed, max_n_segment, height" << std::endl;
    
    for(const auto& item : essential_propagator_codes)
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
    std::cout << "Total number of iterations to compute all propagators: " << total_segments << std::endl;
    std::cout << "------------------------------------" << std::endl;
}

void Molecules::display_sub_propagators() const
{
    // Print sub propagators
    std::vector<std::tuple<std::string, int, int>> sub_deps;
    int total_segments = 0;
    std::cout << "--------- Propagators ---------" << std::endl;
    std::cout << "Key:\n\tsuperposed, max_n_segment, height, deps," << std::endl;
    
    for(const auto& item : essential_propagator_codes)
    {
        total_segments += item.second.max_n_segment;

        std::cout << item.first;
        std::cout << ":\n\t ";
        if (item.first.find('[') == std::string::npos)
            std::cout << "X, ";
        else
            std::cout << "O, ";
        std::cout << item.second.max_n_segment << ", " << item.second.height;

        sub_deps = get_deps_from_key(item.first);
        for(size_t i=0; i<sub_deps.size(); i++)
        {
            std::cout << ", "  << std::get<0>(sub_deps[i]) << ":" << std::get<1>(sub_deps[i]);
        }
        std::cout << std::endl;
    }
    std::cout << "Total number of iterations to compute all propagators: " << total_segments << std::endl;
    std::cout << "------------------------------------" << std::endl;
}

bool ComparePropagatorKey::operator()(const std::string& str1, const std::string& str2)
{
    // First compare heights
    int height_str1 = Molecules::get_height_from_key(str1);
    int height_str2 = Molecules::get_height_from_key(str2);

    if (height_str1 < height_str2)
        return true;
    else if(height_str1 > height_str2)
        return false;

    // second compare their strings
    int mix_length = std::min(str1.length(), str2.length());
    for(int i=0; i<mix_length; i++)
    {
        if (str1[i] == str2[i])
            continue;
        else if (str2[i] == '[')
            return true;
        else if (str1[i] == '[')
            return false;
        else if (str2[i] == ']')
            return true;
        else if (str1[i] == ']')
            return false;
        else if (str2[i] == '(')
            return true;
        else if (str1[i] == '(')
            return false;
        else if (str2[i] == ')')
            return true;
        else if (str1[i] == ')')
            return false;
        else
            return str1[i] < str2[i];
    }
    // Third compare their lengths
    return str1.length() < str2.length();
}