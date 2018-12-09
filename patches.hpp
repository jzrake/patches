#pragma once
#include <array>
#include <map>
#include "ndarray.hpp"




// ============================================================================
namespace patches2d {


    class Database;


    // ========================================================================
    enum class Field
    {
        cell_volume,
        cell_coords,
        vert_coords,
        face_area_i,
        face_area_j,
        conserved,
    };


    // ========================================================================
    enum class MeshLocation
    {
        vert,
        cell,
        face_i,
        face_j,
    };


    // ========================================================================
    struct FieldDescriptor
    {
        FieldDescriptor(int num_fields, MeshLocation location);
        int num_fields;
        MeshLocation location;
    };
}




// ============================================================================
class patches2d::Database
{
public:


    // ========================================================================
    using Header = std::map<Field, FieldDescriptor>;
    using Index = std::tuple<int, int, int, Field>; // i, j, level, which
    using Array = nd::array<double, 3>;


    // ========================================================================
    Database(int ni, int nj, Header header);


    /**
     * Insert a deep copy of the given array into the database at the given
     * patch index. Any existing data at that location is overwritten.
     */
    void insert(Index index, Array data);


    /**
     * Erase any patch data at the given index.
     */
    auto erase(Index index);


    /**
     * Merge data into the database at index, with the given weighting factor.
     * Setting rk_factor=0.0 corresponds to overwriting the existing data.
     *
     * An exception is throws if a patch does not already exist at the given patch
     * index. Use insert to create a new patch.
     */
    void commit(Index index, Array data, double rk_factor=0.0);


    /**
     * Return a deep copy of the data at the patch index, padded with the
     * given number of guard zones at each edge of the array. If no data
     * exists at that index, an exception is thrown.
     *
     *          jl
     *      +--------+
     *      |        |
     *  il  |        |  ir
     *      |        |
     *      +--------+
     *          jr
     *
     */
    Array fetch(Index index, int ngil, int ngir, int ngjl, int ngjr) const;


    /**
     * Same as above, where the number of guard zones to be fetched is the
     * same on each of the patch boundaries.
     */
    Array fetch(Index index, int guard) const;


    /**
     * Return a constant reference to the data at the given patch index. If no
     * data exists at that index, an exception is thrown.
     */
    const Array& at(Index index) const;


    /**
     * Same as above, except discards the index field and uses the given field
     * instead.
     */
    const Array& at(Index index, Field which) const;


    /** Return all patches registered for the given field. */
    std::map<Index, Array> all(Field which) const;

    /** Return an iterator to the beginning of the container of patches. */
    auto begin() const { return patches.begin(); }

    /** Return an iterator to the end of the container of patches. */
    auto end() const { return patches.end(); }

    /** Return the number of patches. */
    std::size_t size() const { return patches.size(); }

    /** Return the number of patches associated with the given field. */
    std::size_t count(Field which) const;

    /** Return the total number of cells associated with the given field. */
    std::size_t num_cells(Field which) const;

    /** Print a description of the patch locations. */
    void print(std::ostream& os) const;

private:
    // ========================================================================
    int num_fields(Index index) const;
    MeshLocation location(Index index) const;
    std::array<int, 3> expected_shape(Index index) const;
    std::array<Index, 4> refine(Index index) const;
    Index coarsen(Index index) const;
    Array check_shape(Array& array, Index index) const;
    Array locate(Index index) const;
    Array quadrant(const nd::array<double, 3>& A, int I, int J) const;
    Array tile(std::array<Index, 4> indexes) const;
    Array prolongation(const nd::array<double, 3>& A) const;
    Array restriction(const nd::array<double, 3>& A) const;

    template <typename IndexContainer>
    bool contains_all(IndexContainer indexes) const
    {
        for (auto index : indexes)
        {
            if (patches.count(index) == 0)
            {
                return false;
            }
        }
        return true;
    }

    // ========================================================================
    int ni;
    int nj;
    Header header;
    std::map<Index, Array> patches;
};




// ============================================================================
namespace patches2d {
    std::string to_string(Database::Index index);
}
