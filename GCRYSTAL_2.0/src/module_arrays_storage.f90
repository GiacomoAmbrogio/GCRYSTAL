Module module_arrays_storage
Use module_arrays

  Implicit None
  Public

  Type(ks_array)   , Public :: base
  Type(ks_array)   , Public :: sk_array
  Type(ks_array)   , Public :: hk_array
  Type(ks_array)   , Public :: fk_array
  Type(ks_array_1D), Public :: evals_array

End Module module_arrays_storage
