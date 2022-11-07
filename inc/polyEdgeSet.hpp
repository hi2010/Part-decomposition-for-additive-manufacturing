/**
 * @file polyEdgeSet.hpp
 * @author your name (you@domain.com)
 * @brief a set to iterate ovet a unique set of edges of a poly data. Currently not used.
 *  How to use: create a UnorederedEdgeSet. insert edges using ?insert?(isertEdgeSorted(ip0, ip1)).
 *  This way it is ensured, that the edge p1 -> p0 does not get inserted if p0 -> p1 already exists
 * @version 0.1
 * @date 2022-03-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef POLY_EDGE_SET_HPP
#define POLY_EDGE_SET_HPP

#include <vtkType.h>
#include <unordered_set>

// id based edge that points to the points in a polydata that belong to one edge
class PolyEdge {
    public:
    void insertEdgeSorted(vtkIdType ip0, vtkIdType ip1){
        if (ip0 < ip1) {
            this->mIdx0 = ip0;
            this->mIdx1 = ip1;
        } else {
            this->mIdx0 = ip1;
            this->mIdx1 = ip0;
        }
    }

    /**
     * @brief checks if the indexes of the edges are the same (mIdx0 / 1 == -> true)
     * 
     * @param e0 edge 0
     * @param e1 edge 1
     * @return true 
     * @return false 
     */
    static bool equals(const PolyEdge& e0, const PolyEdge& e1) {
        if ((e0.mIdx0 == e1.mIdx0) && (e0.mIdx1 == e1.mIdx1)) return true;
        return false;
    }

    /**
     * @brief checks if the indexes of the edges are the same (mIdx0 / 1 == -> true)
     * 
     * @param e the edge to compre this to 
     * @return true 
     * @return false 
     */
    bool equals(const PolyEdge& e) {
        if ((this->mIdx0 == e.mIdx0) && (this->mIdx1 == e.mIdx1)) return true;
        return false;
    }

    vtkIdType getHash(void) {
        // creates a hash, there is some chance that hashes are same for different data (if the data is to big / the idx)
        // -> to big meaning sqrt(vtkIdType)
        // shift the value by the half bit with of the valuetype to the left so that it starts right after the middle.
        auto halflength = (((sizeof(vtkIdType)/sizeof(uint8_t)) >> 1)*8)-1;
        return (this->mIdx1 << (halflength*8)) + this->mIdx0;
    }

    static vtkIdType getHash(const PolyEdge& e) {
        // creates a hash, there is some chance that hashes are same for different data (if the data is to big / the idx)
        // -> to big meaning sqrt(vtkIdType)
        // shift the value by the half bit with of the valuetype to the left so that it starts right after the middle.
        auto halflength = (((sizeof(vtkIdType)/sizeof(uint8_t)) >> 1)*8)-1;
        return (e.mIdx1 << (halflength*8)) + e.mIdx0;
    }

    vtkIdType mIdx0 = 0;
    vtkIdType mIdx1 = 0;
};

// Custom Hash Functor that will compute the hash on the
// passed string objects length
struct EdgeHashByHalfsizeValue {
public:
    vtkIdType operator()(const PolyEdge & e) const {
        return PolyEdge::getHash(e);
    }
};
// Custom comparator that compares the string objects by length
struct EdgeEqualByIdxs {
public:
    bool operator()(const PolyEdge & e0, const PolyEdge & e1) const {
        return PolyEdge::equals(e0, e1);
    }
};

typedef std::unordered_set<PolyEdge, EdgeHashByHalfsizeValue, EdgeEqualByIdxs> UnorderedEdgeSet;

#endif