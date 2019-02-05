#include <ostream>
#include <vector>
#include <map>
#include "patches.hpp"

using namespace patches2d;




// ============================================================================
std::string patches2d::to_string(MeshLocation location)
{
    switch (location)
    {
        case MeshLocation::vert: return "vert";
        case MeshLocation::cell: return "cell";
        case MeshLocation::face_i: return "face_i";
        case MeshLocation::face_j: return "face_j";
    }
}

std::string patches2d::to_string(Field field)
{
    switch (field)
    {
        case Field::cell_volume: return "cell_volume";
        case Field::cell_coords: return "cell_coords";
        case Field::vert_coords: return "vert_coords";
        case Field::face_area_i: return "face_area_i";
        case Field::face_area_j: return "face_area_j";
        case Field::face_velocity_i: return "face_velocity_i";
        case Field::face_velocity_j: return "face_velocity_j";
        case Field::conserved: return "conserved";
        case Field::primitive: return "primitive";
    }
}

std::string patches2d::to_string(Database::Index index)
{
    return to_string(index, to_string(std::get<3>(index)));
}

std::string patches2d::to_string(Database::Index index, std::string field_name)
{
    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto p = std::get<2>(index);
    return std::to_string(p) + "." + std::to_string(i) + "-" + std::to_string(j) + "/" + field_name;
}

MeshLocation patches2d::parse_location(std::string str)
{
    if (str == "vert") return MeshLocation::vert;
    if (str == "cell") return MeshLocation::cell;
    if (str == "face_i") return MeshLocation::face_i;
    if (str == "face_j") return MeshLocation::face_j;
    throw std::invalid_argument("unknown location: " + str);
}

Field patches2d::parse_field(std::string str)
{
    if (str == "cell_volume") return Field::cell_volume;
    if (str == "cell_coords") return Field::cell_coords;
    if (str == "vert_coords") return Field::vert_coords;
    if (str == "face_area_i") return Field::face_area_i;
    if (str == "face_area_j") return Field::face_area_j;
    if (str == "face_velocity_i") return Field::face_velocity_i;
    if (str == "face_velocity_j") return Field::face_velocity_j;
    if (str == "conserved")   return Field::conserved;
    if (str == "primitive")   return Field::primitive;
    throw std::invalid_argument("unknown field: " + str);
}

Database::Index patches2d::parse_index(std::string str)
{
    auto dot   = str.find('.');
    auto dash  = str.find('-');
    auto slash = str.find('/');
    int level  = std::stoi(str.substr(0, dot));
    int i      = std::stoi(str.substr(dot  + 1, dash));
    int j      = std::stoi(str.substr(dash + 1, slash));
    auto field = parse_field(str.substr(slash + 1));
    return std::make_tuple(i, j, level, field);
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

void Database::set_boundary_value(BoundaryValue b)
{
    boundary_value = b;
}

void Database::insert(Index index, Array data)
{
    patches[index].become(check_shape(data, index).copy());
}

auto Database::erase(Index index)
{
    return patches.erase(index);
}

void Database::clear()
{
    patches.clear();
}

void Database::commit(Index index, Array data, double rk_factor)
{
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

Database::Array Database::fetch(Index index, int guard) const
{
    return fetch(index, guard, guard, guard, guard);
}

Database::Array Database::fetch(Index index, int ngil, int ngir, int ngjl, int ngjr) const
{
    if (location(index) != MeshLocation::cell)
    {
        throw std::invalid_argument("can only fetch cell data (for now)");
    }

    auto _     = nd::axis::all();
    auto mi    = ni + ngil + ngir;
    auto mj    = nj + ngjl + ngjr;
    auto shape = std::array<int, 3>{mi, mj, num_fields(index)};
    auto res   = nd::array<double, 3>(shape);

    auto i = std::get<0>(index);
    auto j = std::get<1>(index);
    auto p = std::get<2>(index);
    auto f = std::get<3>(index);
    auto Ri = std::make_tuple(i + 1, j, p, f);
    auto Li = std::make_tuple(i - 1, j, p, f);
    auto Rj = std::make_tuple(i, j + 1, p, f);
    auto Lj = std::make_tuple(i, j - 1, p, f);

    res.select(_|ngil|ni+ngil, _|ngjl|nj+ngjl, _) = patches.at(index);

    // i-left boundary
    // ========================================================================
    if (ngil > 0)
    {
        auto neighbor = locate(Li);
        auto bv = neighbor.empty()
        ? boundary_value(index, PatchBoundary::il, ngil, patches.at(index))
        : neighbor.select(_|ni-ngil|ni, _, _);
        res.select(_|0|ngil, _|ngjl|nj+ngjl, _) = bv;
    }

    // i-right boundary
    // ========================================================================
    if (ngir > 0)
    {
        auto neighbor = locate(Ri);
        auto bv = neighbor.empty()
        ? boundary_value(index, PatchBoundary::ir, ngir, patches.at(index))
        : neighbor.select(_|0|ngir, _, _);
        res.select(_|mi-ngir|mi, _|ngjl|nj+ngjl, _) = bv;
    }

    // j-left boundary
    // ========================================================================
    if (ngjl > 0)
    {
        auto neighbor = locate(Lj);
        auto bv = neighbor.empty()
        ? boundary_value(index, PatchBoundary::jl, ngjl, patches.at(index))
        : neighbor.select(_, _|nj-ngjl|nj, _);
        res.select(_|ngil|ni+ngil, _|0|ngjl, _) = bv;
    }

    // j-right boundary
    // ========================================================================
    if (ngjr > 0)
    {
        auto neighbor = locate(Rj);
        auto bv = neighbor.empty()
        ? boundary_value(index, PatchBoundary::jr, ngjr, patches.at(index))
        : neighbor.select(_, _|0|ngjr, _);
        res.select(_|ngil|ni+ngil, _|mj-ngjr|mj, _) = bv;
    }

    return res;
}

Database::Array Database::assemble(int i0, int i1, int j0, int j1, int level, Field field) const
{
    auto _ = nd::axis::all();
    int mi = 0, mj = 0;

    switch (header.at(field).location)
    {
        case MeshLocation::cell  : mi = (i1 - i0) * ni + 0; mj = (j1 - j0) * nj + 0; break;
        case MeshLocation::vert  : mi = (i1 - i0) * ni + 1; mj = (j1 - j0) * nj + 1; break;
        case MeshLocation::face_i: mi = (i1 - i0) * ni + 1; mj = (j1 - j0) * nj + 0; break;
        case MeshLocation::face_j: mi = (i1 - i0) * ni + 0; mj = (j1 - j0) * nj + 1; break;
    }

    auto res = Array(mi, mj, header.at(field).num_fields);

    for (int i = i0; i < i1; ++i)
    {
        for (int j = j0; j < j1; ++j)
        {
            int di = 0, dj = 0;

            switch (header.at(field).location)
            {
                case MeshLocation::cell  : di = (i == i1 - 1) * 0; dj = (j == j1 - 1) * 0; break;
                case MeshLocation::vert  : di = (i == i1 - 1) * 1; dj = (j == j1 - 1) * 1; break;
                case MeshLocation::face_i: di = (i == i1 - 1) * 1; dj = (j == j1 - 1) * 0; break;
                case MeshLocation::face_j: di = (i == i1 - 1) * 0; dj = (j == j1 - 1) * 1; break;
            }

            const auto& patch = patches.at(std::make_tuple(i, j, level, field));
            res.select(_|i*ni|(i+1)*ni+di, _|j*nj|(j+1)*nj+dj, _) = patch.select(_|0|ni+di, _|0|nj+dj, _);
        }
    }
    return res;
}

const Database::Array& Database::at(Index index) const
{
    return patches.at(index);
}

const Database::Array& Database::at(Index index, Field which) const
{
    std::get<3>(index) = which;
    return patches.at(index);
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
    os << "Database:\n\n";
    os << "block size: " << ni << " " << nj << "\n";
    os << "mesh patches:\n\n";

    for (const auto& patch : patches)
    {
        os << "\t" << to_string(patch.first) << "\n";
    }
    os << "\n";
}

void Database::dump(const Serializer& ser) const
{
    ser.write_header(header);
    ser.write_block_size({ni, nj});

    for (const auto& patch : patches)
    {
        ser.write_array(to_string(patch.first), patch.second);
    }
}

Database Database::load(const Serializer& ser, std::set<Field> fields, std::function<bool()> bailout)
{
    auto header = ser.read_header();
    auto blocks = ser.read_block_size();
    auto database = Database(blocks[0], blocks[1], header);

    for (auto patch : ser.list_patches())
    {
        for (auto field : ser.list_fields(patch))
        {
            if (fields.empty() || fields.count(parse_field (field)))
            {
                auto ind = patch + "/" + field;
                database.insert(parse_index(ind), ser.read_array(ind));

                if (bailout && bailout())
                {
                    return database;
                }
            }
        }
    }
    return database;
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
        case MeshLocation::face_i: return {ni + 1, nj + 0, num_fields(index)};
        case MeshLocation::face_j: return {ni + 0, nj + 1, num_fields(index)};
    }
    throw;
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

    return nd::array<double, 3>();
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
