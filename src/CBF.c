#include "CBF_driver.h"
#include "CBF_parameters.h"

int main(int argc, char* argv[])
{
    // Create parameter data
    UserData data; // instance pointer.
    data = CBF_parameters();

    int retval;
    retval = CBF_driver(data);

    return retval;
}