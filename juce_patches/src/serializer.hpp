#include "JuceHeader.h"
#include <array>
#include <vector>
#include <string>




// ============================================================================
class FileSystemSerializer : public patches2d::Serializer
{
public:
    static bool looksLikeDatabase (juce::File path);
    FileSystemSerializer (juce::File chkpt);
    std::vector<std::string> list_fields (std::string patch_index) const override;
    std::vector<std::string> list_patches() const override;
    nd::array<double, 3> read_array (std::string path) const override;
    std::array<int, 2> read_block_size() const override;
    patches2d::Database::Header read_header() const override;
    void write_array (std::string path, const nd::array<double, 3>& patch) const override;
    void write_header (patches2d::Database::Header header) const override;
    void write_block_size (std::array<int, 2> block_size) const override;

private:
    juce::File chkpt;
};
