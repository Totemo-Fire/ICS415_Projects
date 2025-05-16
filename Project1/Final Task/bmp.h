#include <iostream>
#include <fstream>
#include <vector>

#pragma pack(push, 1)
struct BMPHeader {
    uint16_t fileType{0x4D42}; // "BM"
    uint32_t fileSize{0};
    uint16_t reserved1{0};
    uint16_t reserved2{0};
    uint32_t offsetData{0};
};

struct BMPInfoHeader {
    uint32_t size{0};
    int32_t width{0};
    int32_t height{0};
    uint16_t planes{1};
    uint16_t bitCount{24};
    uint32_t compression{0};
    uint32_t sizeImage{0};
    int32_t xPixelsPerMeter{0};
    int32_t yPixelsPerMeter{0};
    uint32_t colorsUsed{0};
    uint32_t colorsImportant{0};
};

struct BMPColorHeader {
    uint32_t redMask{0x00ff0000};
    uint32_t greenMask{0x0000ff00};
    uint32_t blueMask{0x000000ff};
    uint32_t alphaMask{0xff000000};
    uint32_t colorSpaceType{0x73524742}; // "sRGB"
    uint32_t unused[16]{0};
};
#pragma pack(pop)

void writeBMP(const char* filename, int width, int height, const std::vector<uint8_t>& image) {
    BMPHeader header;
    BMPInfoHeader infoHeader;
    header.fileSize = sizeof(BMPHeader) + sizeof(BMPInfoHeader) + image.size();
    header.offsetData = sizeof(BMPHeader) + sizeof(BMPInfoHeader);
    infoHeader.size = sizeof(BMPInfoHeader);
    infoHeader.width = width;
    infoHeader.height = height;
    infoHeader.sizeImage = image.size();

    std::ofstream file(filename, std::ios::binary);
    if (file) {
        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
        file.write(reinterpret_cast<const char*>(&infoHeader), sizeof(infoHeader));
        file.write(reinterpret_cast<const char*>(image.data()), image.size());
		file.close();
    }
}