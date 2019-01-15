// It contains LUT for lfsr8 counter

#include "clicpix2_frameDecoder.hpp"

using namespace caribou;

const uint8_t clicpix2_frameDecoder::lfsr8_lut[255] = {
    0,   1,   198, 2,   104, 199, 238, 3,   193, 105, 113, 200, 141, 239, 14,  4,   194, 181, 106, 185, 114, 228, 201, 120,
    142, 147, 240, 69,  15,  36,  5,   47,  195, 110, 182, 66,  107, 84,  186, 133, 115, 87,  229, 212, 202, 159, 121, 21,
    189, 143, 208, 148, 92,  241, 56,  70,  136, 16,  152, 37,  98,  6,   221, 48,  196, 236, 111, 12,  183, 118, 67,  45,
    108, 82,  85,  157, 187, 90,  134, 96,  234, 116, 80,  88,  232, 230, 174, 213, 176, 203, 169, 160, 215, 122, 245, 22,
    252, 190, 178, 144, 63,  209, 205, 149, 42,  93,  171, 242, 60,  57,  162, 71,  137, 32,  17,  217, 153, 165, 38,  124,
    99,  128, 7,   247, 222, 74,  49,  24,  254, 197, 103, 237, 192, 112, 140, 13,  180, 184, 227, 119, 146, 68,  35,  46,
    109, 65,  83,  132, 86,  211, 158, 20,  188, 207, 91,  55,  135, 151, 97,  220, 235, 11,  117, 44,  81,  156, 89,  95,
    233, 79,  231, 173, 175, 168, 214, 244, 251, 177, 62,  204, 41,  170, 59,  161, 31,  216, 164, 123, 127, 246, 73,  23,
    253, 102, 191, 139, 179, 226, 145, 34,  64,  131, 210, 19,  206, 54,  150, 219, 10,  43,  155, 94,  78,  172, 167, 243,
    250, 61,  40,  58,  30,  163, 126, 72,  101, 138, 225, 33,  130, 18,  53,  218, 9,   154, 77,  166, 249, 39,  29,  125,
    100, 224, 129, 52,  8,   76,  248, 28,  223, 51,  75,  27,  50,  26,  25};