/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * generate_headers.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#include <fstream>
#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <vector>

bool read(const std::string& name, std::string& s) {
    std::ifstream in(name.c_str());
    
    if (!in) {
        std::cerr << "Could not open " << name << " for reading." << std::endl;
        return false;
    }
    
    std::ostringstream buffer;
    buffer << in.rdbuf(); 

    s = buffer.str();
    
    return true;
}

unsigned replace(const std::string& placeholder, const std::string& replacement, std::string& s) {
    std::string::size_type off = 0;
    unsigned count = 0;
    while ((off = s.find(placeholder, off)) != std::string::npos) {
        s.replace(off, placeholder.size(), replacement);
        ++count;
        ++off;
    }
    
    if (count > 0)
        std::cout << "Replaced \"" << placeholder << "\" with \"" << replacement << "\" " << count << " times." << std::endl;
    return count;
}

bool process(const unsigned dim, const std::string& inName, const std::string& outName) {
    
    std::string s;
    if (!read(inName, s))
        return false;
    
    // Replace all the dimensional placeholders with the actual dimension.
    {
        std::ostringstream replacement;
        replacement << dim;
        replace("$D", replacement.str(), s);
    }
    
    // Replace the vector constructor arguments with the positional arguments.
    {
        std::ostringstream replacement;
        for (unsigned i = 0; i < dim; ++i) 
            replacement << "const T v" << i << (i < dim - 1 ? ", " : "");
            
        replace("$VECTOR_CONSTRUCTOR_ARGUMENTS", replacement.str(), s);
    }
 
    // Replace the vector constructor body.
    {
        std::ostringstream replacement;
        for (unsigned i = 0; i < dim; ++i) 
            replacement << "m_data[" << i << "] = v" << i << "; ";
        
        replace("$VECTOR_CONSTRUCTOR_BODY", replacement.str(), s);
    }
    
    // Replace the matrix constructor arguments with the positional arguments.
    {
        std::ostringstream replacement;
        for (unsigned j = 0; j < dim; ++j)
            for (unsigned i = 0; i < dim; ++i) 
                replacement << "const T v" << i << j << (!(i == dim-1 && j == dim-1) ? ", " : "");
        
        replace("$MATRIX_CONSTRUCTOR_ARGUMENTS", replacement.str(), s);
    }
    
    // Replace the matrix constructor body.
    {
        std::ostringstream replacement;
        for (unsigned j = 0; j < dim; ++j) 
            for (unsigned i = 0; i < dim; ++i) 
                replacement << "m_data[" << (j * dim + i) << "] = v" << i << j << "; ";
        
        replace("$MATRIX_CONSTRUCTOR_BODY", replacement.str(), s);
    }
    
    
    std::ofstream out(outName.c_str());
    out << "// Warning: this file was generated from " << inName << ".  Do not edit if you value your changes!\n\n";
    out << s;
    
    return true;
}

struct FileInfo {
    unsigned dim;
    std::string file, templatePath, outputPath;
};

int main(int argc, char* argv[]) {
    
    std::string inPath = ".", headerPath = ".", bodyPath = ".";
    if (argc > 1) 
        inPath = argv[1];

    if (argc > 2) {
        headerPath = argv[2];
        bodyPath = headerPath;
    }
    
    if (argc > 3) 
        bodyPath = argv[3];
    
    std::cout << "Path to template files is \"" << inPath << "\"" << std::endl;
    std::cout << "Path to output header files is \"" << headerPath << "\"" << std::endl;
    std::cout << "Path to output body files is \"" << bodyPath << "\"" << std::endl;

    
    std::vector<FileInfo> fileInfo;
    for (unsigned i = 2; i <= 4; ++i) {
        FileInfo info;
        
        info.dim = i;
        
        std::ostringstream dim;
        dim << i;
        
        info.file = "Vec" + dim.str() + "Fwd.h";
        info.templatePath = inPath + "/VecFwdTemplate.h";
        info.outputPath = headerPath + "/" + info.file;
        fileInfo.push_back(info);

        info.file = "Vec" + dim.str() + ".h";
        info.templatePath = inPath + "/VecTemplate.h";
        info.outputPath = headerPath + "/" + info.file;
        fileInfo.push_back(info);

        info.file = "Vec" + dim.str() + ".cpp";
        info.templatePath = inPath + "/VecTemplate.cpp";
        info.outputPath = bodyPath + "/" + info.file;
        fileInfo.push_back(info);
        
        info.file = "Mat" + dim.str() + "Fwd.h";
        info.templatePath = inPath + "/MatFwdTemplate.h";
        info.outputPath = headerPath + "/" + info.file;
        fileInfo.push_back(info);
        
        info.file = "Mat" + dim.str() + ".h";
        info.templatePath = inPath + "/MatTemplate.h";
        info.outputPath = headerPath + "/" + info.file;
        fileInfo.push_back(info);
        
        info.file = "Mat" + dim.str() + ".cpp";
        info.templatePath = inPath + "/MatTemplate.cpp";
        info.outputPath = bodyPath + "/" + info.file;
        fileInfo.push_back(info);
    }

    for (unsigned i = 0; i < fileInfo.size(); ++i) {
        std::cout << "Creating " << fileInfo[i].file << std::endl;

        if (!process(fileInfo[i].dim, fileInfo[i].templatePath, fileInfo[i].outputPath)) {
            std::cerr << "Failed!" << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}