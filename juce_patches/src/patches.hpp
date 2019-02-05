#pragma once
#include <array>
#include <map>
#include <set>




// ============================================================================
namespace patches2d {


    class Database;
    class Serializer;


    // ========================================================================
    enum class Field
    {
        cell_volume,
        cell_coords,
        vert_coords,
        face_area_i,
        face_area_j,
        face_velocity_i,
        face_velocity_j,
        conserved,
        primitive,
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
    enum class PatchBoundary
    {
        il, ir, jl, jr,
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


    /**
     * A callback to be invoked when a target patch's guard zone region cannot
     * be calculated from its neighbor patches. The callback receives the
     * index of the target patch (the one whose boundary values are required),
     * a flag indicating which edge of that patch is needed, the depth of the
     * guard zone region, and the data currenly in the target patch. The
     * callback must return an array whose shape matches the patch data, but
     * having the number of guard zones (depth) in the off-bounds axis. For
     * example, if edge = PatchBoundary::il and depth = 2, then the callback
     * must return an array with shape [2, data.shape(1), data.shape(2)].
     * 
     * Use the set_boundary_value method to set the callback. If no callback
     * has been supplied and a call to fetch would require it, then an
     * exception is raised.
     */
    using BoundaryValue = std::function<Array(
        
        Index index,        /**< Index of the target patch */

        PatchBoundary edge, /**< Which edge of that patch to return guard zones for */

        int depth,          /**< The number of guard zones needed */

        const Array& patch  /**< The data in the target patch */

        )>;


    // ========================================================================
    Database(int ni, int nj, Header header);


    /**
     * Set the callback to be invoked when a target patch's guard zone region
     * cannot be found in neighboring patches.
     */
    void set_boundary_value(BoundaryValue);


    /**
     * Insert a deep copy of the given array into the database at the given
     * patch index. Any existing data at that location is overwritten.
     */
    void insert(Index index, Array data);


    /**
     * Erase any patch data at the given index.
     */
    auto erase(Index index);


    /** Clear all of the stored patches from the database. */
    void clear();


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
     * Return an array spanning a rectangular range of blocks at a fixed
     * level. All of the enclosed patches must exist in the database. The
     * upper indexes are non-inclusive. Note that vertex and face data
     * formally has redundancies at the patch boundaries (the right faces of
     * patch i and the left faces of patch i + 1 are the same physical
     * things). This function does not make any attempt to reconcile
     * differences between the data at the redundant locations (it replaces a
     * patch's data values on the right if they can be read from the right
     * neighbor patch), so it's up to you to make sure it's consistent.
     */
    Array assemble(int i0, int i1, int j0, int j1, int level, Field field) const;


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


    /** Write the database using the given serialization scheme. */
    void dump(const Serializer&) const;


    /**
     * Load a database using the given serialization scheme. If the fields
     * argument is empty, then all fields are loaded. Otherwise, only those
     * fields are loaded and returned.
     */
    static Database load(const Serializer&, std::set<Field> fields={}, std::function<bool()> bailout=nullptr);


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
    int ni = 0;
    int nj = 0;
    Header header;
    std::map<Index, Array> patches;
    BoundaryValue boundary_value = nullptr;
};




// ============================================================================
namespace patches2d {
    std::string to_string(MeshLocation location);
    std::string to_string(Field field);
    std::string to_string(Database::Index index);
    std::string to_string(Database::Index index, std::string field_name);
    MeshLocation    parse_location(std::string str);
    Field           parse_field(std::string str);
    Database::Index parse_index(std::string str);
}




// ============================================================================
class patches2d::Serializer
{
public:
    /**
     * Destructor.
     */
    virtual ~Serializer() {}

    /**
     * This method must return a list of the fields in a given patch. If this
     * reader is filesystem-based, it probably returns a list of the file
     * names in the directory for that patch. If it is HDF5, it might return
     * the names of the datasets in the group for that patch..
     */
    virtual std::vector<std::string> list_fields(std::string patch_index) const = 0;

    /**
     * This method must return a list of the patches in the database.
     */
    virtual std::vector<std::string> list_patches() const = 0;

    /**
     * This method must return an array read from the given location.
     */
    virtual nd::array<double, 3> read_array(std::string path) const = 0;

    /**
     * This method must return the block size (ni, nj) for the database at the
     * given path.
     */
    virtual std::array<int, 2> read_block_size() const = 0;

    /**
     * This method must return the header for the database at the given path.
     */
    virtual Database::Header read_header() const = 0;

    /**
     * This method must write an array of patch data to the given location.
     */
    virtual void write_array(std::string path, const nd::array<double, 3>& patch) const = 0;

    /**
     * This method must write a header to the given location.
     */
    virtual void write_header(Database::Header header) const = 0;

    /**
     * This method must write a block size (ni, nj) to the given location.
     */
    virtual void write_block_size(std::array<int, 2> block_size) const = 0;
};
