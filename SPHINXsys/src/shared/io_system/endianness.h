/**
 * @file 	endianness.cpp
 * @brief 	Classes for system endianness.
 * @author	Constantinos Menelaou
 */
#pragma once
#include <fstream>

namespace SPH
{
    /**
     * @class Endianness
     * @brief Enum class for endianness of the system.
     */
    enum class Endianness
    {
        little = 0,
        big = 1
    };

    /**
     * @class Endianness
     * @brief Enum class for endianness of the system.
     */
    class Endian
    {
    public:
        /**
         * @brief Get the system endianness.
         * @return Endianness::big or Endianness::little.
         */
        static Endianness getSystemEndianness();

        /**
         * @brief Check if the system is big endian.
         * @return True if the system is big endian.
         */
        static bool isBigEndian();
        
        /**
         * @brief Check if the system is little endian.
         * @return True if the system is little endian.
         */
        static bool isLittleEndian();
        
        /**
         * @brief Write data to file in reverse endianness.
         * @param file The file to write to.
         * @param data The data to write.
         * @param bytes_count The number of bytes to write.
         * @param type_size The size of the data type.
         */
        static void writeDataReverseEndianness(
            std::ofstream &file, const void *data, size_t bytes_count, size_t type_size);
    };
} // namespace SPH
