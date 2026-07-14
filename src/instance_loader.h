#pragma once

#include <filesystem>
#include "model.h"

Model load_instance_from_flat(const std::filesystem::path& file_path);
