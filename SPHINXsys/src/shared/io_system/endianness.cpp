/**
 * @file 	endianness.cpp
 * @brief 	Definitions for Endian class.
 * @author	Constantinos Menelaou
 */
#include "endianness.h"

namespace SPH
{
    //=============================================================================================//
    Endianness Endian::getSystemEndianness()
    {
        const int value  = 0x01;
        const void *address = static_cast<const void *>(&value);
        const unsigned char *least_significant_address = static_cast<const unsigned char *>(address);
        return (*least_significant_address == 0x01) ? Endianness::little : Endianness::big;
    }
    //=============================================================================================//
    bool Endian::isBigEndian()
    {
        return (Endian::getSystemEndianness() == Endianness::big);
    }
    //=============================================================================================//
    bool Endian::isLittleEndian()
    {
        return (Endian::getSystemEndianness() == Endianness::little);
    }
    //=============================================================================================//
    void Endian::writeDataReverseEndianness(std::ofstream &file, const void *data, size_t bytes_count, size_t type_size)
    {
        const char *data_char = reinterpret_cast<const char *>(data);
        for (int i = 0; i < bytes_count; i += type_size) {
            for (int j = type_size-1; j >=0; j--) {
                file.write((data_char+i+j), 1);
            }
        }
    }
    //=============================================================================================//    
} // namespace SPH

