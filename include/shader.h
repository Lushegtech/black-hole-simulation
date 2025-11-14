#ifndef SHADER_H
#define SHADER_H

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

/**
 * Shader utility class for loading and compiling GLSL shaders
 */
class Shader {
  public:
    unsigned int ID; // Program ID

    /**
     * Constructor: load and compile shaders from file paths
     */
    Shader(const char* vertexPath, const char* fragmentPath);

    /**
     * Activate the shader
     */
    void use() const;

    /**
     * Utility uniform functions
     */
    void setBool(const std::string& name, bool value) const;
    void setInt(const std::string& name, int value) const;
    void setFloat(const std::string& name, float value) const;
    void setVec2(const std::string& name, float x, float y) const;
    void setVec3(const std::string& name, float x, float y, float z) const;
    void setVec4(const std::string& name, float x, float y, float z, float w) const;

  private:
    /**
     * Check compilation/linking errors
     */
    void checkCompileErrors(unsigned int shader, std::string type);
};

#endif // SHADER_H
