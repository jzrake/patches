#include <ostream>
#include "patches.hpp"

using namespace patches2d;




// ============================================================================
std::string patches2d::to_string(Database::Index index)
{
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto p = std::get<2>(index);
    auto f = std::string();

    switch (std::get<3>(index))
    {
        case Field::conserved: f = "conserved"; break;
        case Field::cell_coords: f = "cell_coords"; break;
        case Field::vert_coords: f = "vert_coords"; break;
    }
    return std::to_string(p) + ":" + std::to_string(i) + "-" + std::to_string(j) + "/" + f;
}




// ============================================================================
FieldDescriptor::FieldDescriptor(int num_fields, MeshLocation location)
: num_fields(num_fields)
, location(location)
{
}




// ============================================================================
Database::Database(int ni, int nj, Header header)
: ni(ni)
, nj(nj)
, header(header)
{
}

void Database::insert(Index index, Array data)
{
    patches[index].become(check_shape(data, index).copy());
}

auto Database::erase(Index index)
{
    return patches.erase(index);
}

void Database::commit(Index index, Array data, double rk_factor)
{
    if (location(index) != MeshLocation::cell)
    {
        throw std::invalid_argument("Can only commit cell data (for now)");
    }

    auto target = patches.at(index);

    if (rk_factor == 0.0)
    {
        target = data;
    }
    else
    {
        target = data * (1 - rk_factor) + target * rk_factor;
    }
}

Database::Array Database::checkout(Index index, int guard) const
{
    if (location(index) != MeshLocation::cell)
    {
        throw std::invalid_argument("Can only checkout cell data (for now)");
    }

    auto _     = nd::axis::all();
    auto ng    = guard;
    auto shape = std::array<int, 3>{ni + 2 * ng, nj + 2 * ng, 5};
    auto res   = nd::array<double, 3>(shape);

    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto p = std::get<2>(index);
    auto f = std::get<3>(index);

    auto Ri = std::make_tuple(i + 1, j, p, f);
    auto Li = std::make_tuple(i - 1, j, p, f);
    auto Rj = std::make_tuple(i, j + 1, p, f);
    auto Lj = std::make_tuple(i, j - 1, p, f);

    res.select(_|ng|ni+ng, _|ng|nj+ng, _) = patches.at(index);

    res.select(_|ni+ng|ni+2*ng, _|ng|nj+ng, _) = locate(Ri).select(_|0|ng, _, _);
    res.select(_|ng|ni+ng, _|nj+ng|nj+2*ng, _) = locate(Rj).select(_, _|0|ng, _);

    res.select(_|0|ng, _|ng|nj+ng, _) = locate(Li).select(_|ni-ng|ni, _, _);
    res.select(_|ng|ni+ng, _|0|ng, _) = locate(Lj).select(_, _|nj-ng|nj, _);

    return res;
}

std::map<Database::Index, Database::Array> Database::all(Field which) const
{
    auto res = std::map<Index, Array>();

    for (const auto& patch : patches)
    {
        if (std::get<3>(patch.first) == which)
        {
            res.insert(patch);
        }
    }
    return res;
}

std::size_t Database::count(Field which) const
{
    std::size_t n = 0;

    for (const auto& patch : patches)
    {
        if (std::get<3>(patch.first) == which)
        {
            ++n;
        }
    }
    return n;
}

std::size_t Database::num_cells(Field which) const
{
    return count(which) * ni * nj;
}

void Database::print(std::ostream& os) const
{
    os << std::string(52, '=') << "\n";
    os << "Mesh patches:\n\n";

    for (const auto& patch : patches)
    {
        os << to_string(patch.first) << "\n";
    }
    os << "\n";
}




// ========================================================================
int Database::num_fields(Index index) const
{
    return header.at(std::get<3>(index)).num_fields;
}

MeshLocation Database::location(Index index) const
{
    return header.at(std::get<3>(index)).location;
}

Database::Array Database::check_shape(Array& array, Index index) const
{
    if (array.shape() != expected_shape(index))
    {
        throw std::invalid_argument("input patch data has the wrong shape");
    }
    return array;
}

Database::Index Database::coarsen(Index index) const
{
    std::get<0>(index) /= 2;
    std::get<1>(index) /= 2;
    std::get<2>(index) -= 1;
    return index;
}

std::array<Database::Index, 4> Database::refine(Index index) const
{
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto p = std::get<2>(index);
    auto f = std::get<3>(index);

    return {
        std::make_tuple(i * 2 + 0, j * 2 + 0, p + 1, f),
        std::make_tuple(i * 2 + 0, j * 2 + 1, p + 1, f),
        std::make_tuple(i * 2 + 1, j * 2 + 0, p + 1, f),
        std::make_tuple(i * 2 + 1, j * 2 + 1, p + 1, f),
    };
}

std::array<int, 3> Database::expected_shape(Index index) const
{
    switch (location(index))
    {
        case MeshLocation::cell: return {ni + 0, nj + 0, num_fields(index)};
        case MeshLocation::vert: return {ni + 1, nj + 1, num_fields(index)};
    }
}

nd::array<double, 3> Database::locate(Index index) const
{
    if (patches.count(index))
    {
        return patches.at(index);
    }

    if (patches.count(coarsen(index)))
    {
        auto i = std::get<0>(index);
        auto j = std::get<1>(index);
        return prolongation(quadrant(patches.at(coarsen(index)), i % 2, j % 2));
    }

    if (contains_all(refine(index)))
    {
        return restriction(tile(refine(index)));
    }

    // Return a value based on some physical boundary conditions

    auto _ = nd::axis::all();
    auto res = nd::array<double, 3>(ni, nj, num_fields(index));

    res.select(_, _, 0) = 0.1;
    res.select(_, _, 1) = 0.0;
    res.select(_, _, 2) = 0.0;
    res.select(_, _, 3) = 0.0;
    res.select(_, _, 4) = 0.125;

    return res;
}

nd::array<double, 3> Database::quadrant(const nd::array<double, 3>& A, int I, int J) const
{
    auto _ = nd::axis::all();

    if (I == 0 && J == 0) return A.select(_|0 |ni/2, _|0 |nj/2, _);
    if (I == 0 && J == 1) return A.select(_|0 |ni/2, _|nj/2|nj, _);
    if (I == 1 && J == 0) return A.select(_|ni/2|ni, _|0 |nj/2, _);
    if (I == 1 && J == 1) return A.select(_|ni/2|ni, _|nj/2|nj, _);

    throw std::invalid_argument("quadrant: I and J must be 0 or 1");
}

nd::array<double, 3> Database::tile(std::array<Index, 4> indexes) const
{
    auto _ = nd::axis::all();
    nd::array<double, 3> res(ni * 2, nj * 2, num_fields(indexes[0]));

    res.select(_|0 |ni*1, _|0 |nj*1, _) = patches.at(indexes[0]);
    res.select(_|0 |ni*1, _|nj|nj*2, _) = patches.at(indexes[1]);
    res.select(_|ni|ni*2, _|0 |nj*1, _) = patches.at(indexes[2]);
    res.select(_|ni|ni*2, _|nj|nj*2, _) = patches.at(indexes[3]);

    return res;
}

nd::array<double, 3> Database::prolongation(const nd::array<double, 3>& A) const
{
    auto _ = nd::axis::all();
    nd::array<double, 3> res(ni, nj, A.shape(2));

    res.select(_|0|ni|2, _|0|nj|2, _) = A;
    res.select(_|0|ni|2, _|1|nj|2, _) = A;
    res.select(_|1|ni|2, _|0|nj|2, _) = A;
    res.select(_|1|ni|2, _|1|nj|2, _) = A;

    return res;
}

nd::array<double, 3> Database::restriction(const nd::array<double, 3>& A) const
{
    auto _ = nd::axis::all();
    auto mi = A.shape(0);
    auto mj = A.shape(1);

    auto B = std::array<nd::array<double, 3>, 4>
    {
        A.select(_|0|mi|2, _|0|mj|2, _),
        A.select(_|0|mi|2, _|1|mj|2, _),
        A.select(_|1|mi|2, _|0|mj|2, _),
        A.select(_|1|mi|2, _|1|mj|2, _),
    };
    return (B[0] + B[1] + B[2] + B[3]) * 0.25;
}
