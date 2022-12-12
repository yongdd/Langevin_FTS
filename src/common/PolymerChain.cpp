#include <iostream>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <stack>

#include "PolymerChain.h"
#include "Exception.h"

//----------------- Constructor ----------------------------
PolymerChain::PolymerChain(
    double ds, std::map<std::string, double> bond_lengths, 
    double volume_fraction, 
    std::vector<std::string> block_species,
    std::vector<double> contour_lengths,
    std::vector<int> v, std::vector<int> u,
    std::map<int, int> v_to_grafting_index)
{
    // check block size
    if( block_species.size() != contour_lengths.size())
        throw_with_line_number("The sizes of block_species (" + std::to_string(block_species.size()) + 
            ") and contour_lengths (" +std::to_string(contour_lengths.size()) + ") must be consistent.");

    if( block_species.size() != v.size())
        throw_with_line_number("The sizes of block_species (" + std::to_string(block_species.size()) + 
            ") and edges v (" +std::to_string(v.size()) + ") must be consistent.");

    if( block_species.size() != u.size())
        throw_with_line_number("The sizes of block_species (" + std::to_string(block_species.size()) + 
            ") and edges u (" +std::to_string(v.size()) + ") must be consistent.");

    // check block lengths, segments, types
    for(int i=0; i<contour_lengths.size(); i++)
    {
        if( contour_lengths[i] <= 0)
            throw_with_line_number("contour_lengths[" + std::to_string(i) + "] (" +std::to_string(contour_lengths[i]) + ") must be a positive number.");
        if( std::abs(std::lround(contour_lengths[i]/ds)-contour_lengths[i]/ds) > 1.e-6)
            throw_with_line_number("contour_lengths[" + std::to_string(i) + "]/ds (" + std::to_string(contour_lengths[i]) + "/" + std::to_string(ds) + ") is not an integer.");
        if( bond_lengths.count(block_species[i]) == 0 )
            throw_with_line_number("block_species[" + std::to_string(i) + "] (\"" + block_species[i] + "\") is not in bond_lengths.");
    }

    //  
    if( v_to_grafting_index.size() > 0)
        throw_with_line_number("Currently, \'v_to_grafting_index\' is not supported.");

    // save variables
    try
    {
        this->volume_fraction = volume_fraction;
        for(int i=0; i<contour_lengths.size(); i++){
            blocks.push_back({
                block_species[i],                         // species
                (int) std::lround(contour_lengths[i]/ds),   // n_segment
                contour_lengths[i],                         // contour_length
                v[i], u[i]});                               // nodes
        }
    }
    catch(std::exception& exc)
    {
        throw_without_line_number(exc.what());
    }

    // compute alpha, sum of relative contour lengths
    double alpha{0.0};
    for(int i=0; i<blocks.size(); i++){
        alpha += blocks[i].contour_length;
    }
    this->alpha = alpha;

    // construct adjacent_nodes
    for(int i=0; i<contour_lengths.size(); i++)
    {
        adjacent_nodes[v[i]].push_back(u[i]);
        adjacent_nodes[u[i]].push_back(v[i]);
        // v and u must be a non-negative integer for depth first search
        if( v[i] < 0 )
            throw_with_line_number("v[" + std::to_string(i) + "] (" + std::to_string(v[i]) + ") must be a non-negative integer.");
        if( u[i] < 0 )
            throw_with_line_number("u[" + std::to_string(i) + "] (" + std::to_string(u[i]) + ") must be a non-negative integer.");
        if( v[i] == u[i] )
            throw_with_line_number("v[" + std::to_string(i) + "] and u[" + std::to_string(i) + "] must be different integers.");
    }

    // // print adjacent_nodes
    // for(const auto& node : adjacent_nodes){
    //     std::cout << node.first << ": [";
    //     for(int i=0; i<node.second.size()-1; i++)
    //         std::cout << node.second[i] << ",";
    //     std::cout << node.second[node.second.size()-1] << "]" << std::endl;
    // }

    // detect a cycle and isolated nodes in the block copolymer graph using depth first search
    std::map<int, bool> is_visited;
    for (int i = 0; i < contour_lengths.size(); i++)
        is_visited[v[i]] = false;

    std::stack<std::pair<int,int>> connected_nodes;
    connected_nodes.push(std::make_pair(v[0],-1));
    while (!connected_nodes.empty())
    {
        //std::cout << "connected_nodes" << connected_nodes.top() << std::endl;

        // pop item and visit
        int cur = connected_nodes.top().first;
        int parent = connected_nodes.top().second;
        is_visited[cur] = true;
        connected_nodes.pop();

        // add adjacent_nodes at stack
        auto nodes = adjacent_nodes[cur];
        for(int i=0; i<nodes.size();i++)
        {
            if (is_visited[nodes[i]] && nodes[i] != parent)
            {
                throw_with_line_number("A cycle is detected, which contains nodes " 
                    + std::to_string(nodes[i]) + " and " + std::to_string(parent)
                    + ". Only acyclic block copolymer is allowed.");
            }
            else if(! is_visited[nodes[i]])
            {
                connected_nodes.push(std::make_pair(nodes[i], cur));
            }
        }
    }
    for (int i=0; i<contour_lengths.size(); i++)
    {
        if (!is_visited[v[i]])
            throw_with_line_number("There are disconnected nodes. Please check node number: " + std::to_string(v[i]) + ".");
    }

    // construct edge nodes
    for (int i=0; i<contour_lengths.size(); i++)
    {
        if (edge_to_array.count(std::make_pair(v[i], u[i])) > 0)
        {
            throw_with_line_number("There are duplicated edges. Please check the edge between ("
                + std::to_string(v[i]) + ", " + std::to_string(u[i]) + ").");
        }
        else
        {
            edge_to_array[std::make_pair(v[i],u[i])] = i;
            edge_to_array[std::make_pair(u[i],v[i])] = i;
        }
    }
}
int PolymerChain::get_n_blocks()
{
    return blocks.size();
}
int PolymerChain::get_n_segment(int idx)
{
    return blocks[idx].n_segment;
}
int PolymerChain::get_n_segment_total()
{
    int total_n{0};
    for(int i=0; i<blocks.size(); i++)
        total_n += blocks[i].n_segment;
    return total_n;
}
double PolymerChain::get_alpha()
{
    return alpha;
}
double PolymerChain::get_volume_fraction()
{
    return volume_fraction;
}
int PolymerChain::get_array_idx(int v, int u)
{
    if (edge_to_array.count(std::make_pair(v, u)) == 0)
        throw_with_line_number("There is no such edge (" + std::to_string(v) + ", " + std::to_string(u) + ").");
    return edge_to_array[std::make_pair(v, u)];
}
struct PolymerChainBlock& PolymerChain::get_block(int v, int u)
{
    if (edge_to_array.count(std::make_pair(v, u)) == 0)
        throw_with_line_number("There is no such edge (" + std::to_string(v) + ", " + std::to_string(u) + ").");
    return blocks[edge_to_array[std::make_pair(v, u)]];
}
std::vector<PolymerChainBlock>& PolymerChain::get_blocks()
{
    return blocks;
}
std::map<int, std::vector<int>>& PolymerChain::get_adjacent_nodes()
{
    return adjacent_nodes;
}
std::map<std::pair<int, int>, int>& PolymerChain::get_edge_to_array()
{
    return edge_to_array;
}
void PolymerChain::set_edge_to_deps(int v, int u, std::string deps)
{
    edge_to_deps[std::make_pair(v, u)] = deps;
}
std::string PolymerChain::get_dep(int v, int u){
    if (edge_to_deps.count(std::make_pair(v, u)) == 0)
        throw_with_line_number("There is no such block (" + std::to_string(v) + ", " + std::to_string(u) + ").");
    return edge_to_deps[std::make_pair(v,u)];
}