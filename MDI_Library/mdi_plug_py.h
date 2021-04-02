/*! \file
 *
 * \brief Implementation of library-based communication
 */

#ifndef MDI_PLUG_PY
#define MDI_PLUG_PY

#include <Python.h>

/*! \brief Flag whether the Python interpreter has been initialized */
extern int python_interpreter_initialized;

/*! \brief Pointer to the original Python interpreter's dictionary.
 * Only used for Python plugins */
extern PyObject* python_interpreter_dict;

int python_plugin_init( const char* engine_name, const char* engine_path, const char* options, void* engine_comm_ptr );

#endif
