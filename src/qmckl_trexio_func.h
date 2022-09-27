qmckl_exit_code
qmckl_trexio_read(const qmckl_context context,
                  const char* file_name,
                  const int64_t size_max);

#ifdef HAVE_DEVICE_POINTERS
qmckl_exit_code
qmckl_trexio_read_device(const qmckl_context context,
                         const char* file_name,
                         const int64_t size_max,
	                       int device_id);
#endif
