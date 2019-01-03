#include "JuceHeader.h"
#include <fstream>
using namespace patches2d;




// ============================================================================
FileSystemSerializer::FileSystemSerializer (std::string chkpt) : chkpt(chkpt)
{
}

std::vector<std::string> FileSystemSerializer::list_fields (std::string patch_index) const
{
    auto res = std::vector<std::string>();
    auto patch = File (chkpt).getChildFile (patch_index);

    for (auto d : patch.findChildFiles (File::findFiles, false))
    {
        if (d.existsAsFile())
        {
            res.push_back (d.getRelativePathFrom (patch).toStdString());
        }
    }
    return res;
}

std::vector<std::string> FileSystemSerializer::list_patches() const
{
    auto res = std::vector<std::string>();

    for (auto d : File (chkpt).findChildFiles (File::findDirectories, false))
    {
        res.push_back (d.getFileName().toStdString());
    }
    return res;
}

nd::array<double, 3> FileSystemSerializer::read_array (std::string path) const
{
    auto ifs = std::ifstream (File (chkpt).getChildFile (path).getFullPathName().toStdString());
    auto str = std::string (std::istreambuf_iterator<char> (ifs), std::istreambuf_iterator<char>());
    return nd::array<double, 3>::loads (str);
}

std::array<int, 2> FileSystemSerializer::read_block_size() const
{
    auto j = JSON::parse (File (chkpt).getChildFile ("block_size.json"));
    return {j["ni"], j["nj"]};
}

Database::Header FileSystemSerializer::read_header() const
{
    auto j = JSON::parse (File (chkpt).getChildFile ("header.json"));

    if (! j.getDynamicObject())
    {
        throw std::runtime_error ("corrupt database header");
    }
    Database::Header header;

    for (auto field : j.getDynamicObject()->getProperties())
    {
        auto num = int (field.value[0]);
        auto loc = patches2d::parse_location (field.value[1].toString().toStdString());
        auto ind = patches2d::parse_field (field.name.toString().toStdString());
        header.emplace (ind, FieldDescriptor (num, loc));
    }
    return header;
}

void FileSystemSerializer::write_array (std::string path, const nd::array<double, 3>& patch) const
{
    throw std::logic_error ("serializer is read-only");
}

void FileSystemSerializer::write_header (Database::Header header) const
{
    throw std::logic_error ("serializer is read-only");
}

void FileSystemSerializer::write_block_size (std::array<int, 2> block_size) const
{
    throw std::logic_error ("serializer is read-only");
}
