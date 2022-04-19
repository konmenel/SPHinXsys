/**
 * @file 	endianness.inl
 * @brief 	Inline definitions for Endian class.
 * @author	Constantinos Menelaou
 */
#pragma once

//=============================================================================================//
inline SPH::Endianness SPH::Endian::getSystemEndianness()
{
    const int value  = 0x01;
    const void *address = static_cast<const void *>(&value);
    const unsigned char *least_significant_address = static_cast<const unsigned char *>(address);
    return (*least_significant_address == 0x01) ? SPH::Endianness::little : SPH::Endianness::big;
}
//=============================================================================================//
inline bool SPH::Endian::isBigEndian()
{
    return (SPH::Endian::getSystemEndianness() == SPH::Endianness::big);
}
//=============================================================================================//
inline bool SPH::Endian::isLittleEndian()
{
    return (SPH::Endian::getSystemEndianness() == SPH::Endianness::little);
}
//=============================================================================================//
inline void SPH::Endian::writeDataReverseEndianness(std::ofstream &file, const void *data, size_t bytes_count, size_t type_size)
{
    const char *data_char = reinterpret_cast<const char *>(data);
    for (int i = 0; i < bytes_count; i += type_size) {
        for (int j = type_size-1; j >=0; j--) {
            file.write((data_char+i+j), 1);
        }
    }
}
//=============================================================================================//

