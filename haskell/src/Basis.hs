module Basis
  (
  ) where

-- | 3D Spherical Basis.
--
-- The purpose of this module is to allow for easy enumeration of a 3D spherical
-- basis, and filter this basis for certain symmetries.
--
-- The basis elements are of the form
-- \[
--   \phi_{i}(r, \Omega)
--     = \phi_{k_{i}, l_{i}, m_{i}}(r, \Omega)
--     = \trac{1}{r} \varphi_{k_{i}, l_{i}}(r) Y_{l_{i}}^{m_{i}}(\Omega)
-- \]
-- for \(k = 1, 2, ...\), \(l = 0, 1, ...\), \(m = -l, -l+1, ..., l)\.
--
-- [/Basis Indexing/]
-- For given \((l, m)\), we represent the indexed set of basis elements,
-- \(\mathcal{B}^{(l, m)}\), by
-- > bLM = (ks, (l, m)) :: ([Int], (Int, Int))
-- where @ks = [1 .. ]@.
-- Note that the values of @l@ and @m@ should be such that
-- @l >= 0@ and @abs(m) <= l@, for approriate usage.
--
-- For an arbitrary range of \((l, m)\), we represent the set of basis elements,
-- \(\mathcal{B}\), by
-- > b = [(ks, lm) | lm <- lms] :: [([Int], (Int, Int))]
-- >   where
-- >     ks = [1 .. ]
-- >     lms = [concat [(l, m) | m <- [-l .. l]] | l <- [0 .. ]]
-- That is, we construct the total basis set \(\mathcal{B}\) from the union of
-- the indexed basis sets \(\mathcal{B}^{(l, m)}\), ordered lexiographically by
-- increasing \(l\), and by increasing \(m\) for equivalent \(l\).
-- This ordering is preserved under any filtering of basis elements.
--
-- [/Example Usage/]
--
-- To construct an infinite basis:
-- > b = basis
-- To construct an infinite basis, with \(m\) as a symmetry:
-- > b m = filter (\(ks, (l', m')) -> m' == m) basis
-- > b m = filter ((==) m . projection) basis  -- equivalent
-- To construct an infinite basis, with \(l\) as a symmetry:
-- > b l = filter (\(ks, (l', m')) -> l' == l) basis
-- > b l = filter ((==) l . angular) basis  -- equivalent
-- To construct an infinite basis, with parity \(\p\) as a symmetry:
-- > b p = filter (\(ks, (l', m')) -> (-1)^l == p) basis
-- > b p = filter ((==) p . parity) basis  -- equivalent
--
-- To construct a finite basis, with \(l_{\rm{max}} = 8\) and with \(20 - l\)
-- radial basis functions per \(l\):
-- > b = truncateK (\l -> 20 - l) $ truncateL 8 basis
-- To construct a finite basis, with \(l_{\rm{max}} = 8\) and with \(20 - l\)
-- radial basis functions per \(l\), with parity \(pi\) as a symmetry:
-- > b p = truncateK (\l -> 20 - l)
-- >     $ truncateL 8
-- >     $ filter ((==) p . parity) basis
--
-- Finally, to flatten the set of indexed bases into a single set of basis
-- elements:
-- > b = truncateK (\l -> 10) $ truncateL 8 basis
-- > bes = flatten b
-- Note that at least one of @truncateK@ or @truncateL@ should be used prior to
-- flattening (to avoid interacting with a doubly infinite list).


-- | Total basis.
basis :: [([Int], (Int, Int))]
basis = concat [[([1 .. ], (l, m)) | m <- [-l .. l]] | l <- [0 .. ]]


-- | Value of @l@ for indexed basis \(\mathcal{B}^{(l, m)}\).
angular :: ([Int], (Int, Int)) -> Int
angular (ks, (l, m)) = l

-- | Value of @m@ for indexed basis \(\mathcal{B}^{(l, m)}\).
projection :: ([Int], (Int, Int)) -> Int
projection (ks, (l, m)) = m

-- | Value of @(-1)^l@ (parity) for indexed basis \(\mathcal{B}^{(l, m)}\).
parity :: ([Int], (Int, Int)) -> Int
parity (ks, (l, m)) = if odd l then -1 else 1


-- | Truncate each indexed basis \(\mathcal{B}^{(l, m)}\) to @kPerL l@ elements.
truncateK :: (Int -> Int) -> [([Int], (Int, Int))] -> [([Int], (Int, Int))]
truncateK kPerL b = fmap (\(ks, (l, m)) -> (take (kPerL l) ks, (l, m))) b

-- | Truncate the total basis \(\mathcal{B}\) to contain only indexed bases
-- \(\mathcal{B}^{(l, m)}\) for which @l <= lmax@.
truncateL :: Int -> [([Int], (Int, Int))] -> [([Int], (Int, Int))]
truncateL lmax b = takeWhile (\(ks, (l, m)) -> l <= lmax) b


-- | Flatten set of (l, m) indexed radial bases, into one basis.
flatten :: [([Int], (Int, Int))] -> [(Int, Int, Int)]
flatten b = concat $ fmap annotate b
  where
    annotate :: ([Int], (Int, Int)) -> [(Int, Int, Int)]
    annotate (ks, (l, m)) = (fmap (\k -> (k, l, m)) ks)
