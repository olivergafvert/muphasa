# Muphasa

Muphasa is a software to compute presentations of multi-parameter persistent homology. It was devoloped by Matías R. Bender and Oliver Gäfvert, and Michael Lesnick.

## Instalation

To install Muphasa you need CMake. There are two options:

- If you want to install Muphasa, together with its Python package and interface, you should use the script [/scripts/install_python_interface](https://github.com/olivergafvert/muphasa/blob/main/scripts/install_python_interface) . You might need additional permitions.

- If you want to compile only its console version, you can use the following lines, 
   ```
   mkdir build; cd build
   cmake ../mph
   cmake --build .
   ```

## Examples

- If you installed the Python interface, you can try the notebook file [/python/notebooks/examples.ipynb](https://github.com/olivergafvert/muphasa/blob/main/python/notebooks/examples.ipynb) .

- If you want to test the console version, you can run the following lines
  ```
  ./mph ../mph/examples/example1.txt
  ./mph --firep ../mph/examples/example2.firep
  ```
