void MSprintf(const int f, const char *const h, const MS m)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  if ((h ? dprintf(f, "\n%s: ", h) : 0) < 0)
    perror("dprintf 0");

  const unsigned u = MS2U(m);
  for (unsigned i = 0u, o = (1u << (VSL - 1u)); i < VSL; ++i) {
    if (1 != dprintf(f, "%c", ((u & o) ? '1' : '0')))
      perror("dprintf 1");
    o >>= 1u;
  }

  if (8 != dprintf(f, " (%04X)\n", u))
    perror("dprintf 2");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}

void MDprintf(const int f, const char *const h, const MD m)
{
#ifdef _OPENMP
#pragma omp critical
  {
#endif /* _OPENMP */
  if ((h ? dprintf(f, "\n%s: ", h) : 0) < 0)
    perror("dprintf 0");

  const unsigned u = MD2U(m);
  for (unsigned i = 0u, o = (1u << (VDL - 1u)); i < VDL; ++i) {
    if (1 != dprintf(f, "%c", ((u & o) ? '1' : '0')))
      perror("dprintf 1");
    o >>= 1u;
  }

  if (6 != dprintf(f, " (%02X)\n", u))
    perror("dprintf 2");
#ifdef _OPENMP
  }
#endif /* _OPENMP */
}
