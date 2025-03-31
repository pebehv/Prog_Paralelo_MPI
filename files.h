/*#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>

void read_values_from_file(const char * file, float * data, size_t size) {
    std::ifstream values(file, std::ios::binary);
    values.read(reinterpret_cast<char*>(data), size);
    values.close();
}

void write_values_to_file(const char * file, float * data, size_t size) {
    std::ofstream values(file, std::ios::binary);
    values.write(reinterpret_cast<char*>(data), size);
    values.close();
}
*/

#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

void read_values_from_file(const char *file, float *data, size_t size, int Body) {
    FILE *values = fopen(file, "rb");
    if (file == NULL) {
        perror("Error al abrir el archivo");
        return;
    }
    if (values != NULL) {
        //fread(data, sizeof(float), size, values);
        size_t elements_read = fread(data, Body, size, values);
        if (elements_read != size) {
            perror("Error al leer datos del archivo");
        }

        // Imprime los datos leídos
       /* for (size_t i = 0; i < elements_read; i++) {
            printf("::%f\n", data[i]);
        }*/
        
    } else {
        // Manejo de error al abrir el archivo
        perror("Error al abrir el archivo");
    }
    fclose(values);
}

void write_values_to_file(const char *file, float *data, size_t size) {
    FILE *values = fopen(file, "wb");
    if (values != NULL) {
        fwrite(data, sizeof(float), size, values);
        fclose(values);
    } else {
        // Manejo de error al abrir el archivo
        perror("Error al abrir el archivo");
    }
}