/*--

This file is a part of libsais, a library for linear time suffix array,
longest common prefix array and burrows wheeler transform construction.

   Copyright (c) 2021-2025 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information.

--*/

#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut
)]
extern "C" {
    fn libsais16(
        T: *const uint16_t,
        SA: *mut int32_t,
        n: int32_t,
        fs: int32_t,
        freq: *mut int32_t,
    ) -> int32_t;
    fn libsais16_gsa(
        T: *const uint16_t,
        SA: *mut int32_t,
        n: int32_t,
        fs: int32_t,
        freq: *mut int32_t,
    ) -> int32_t;
    fn libsais16_int(
        T: *mut int32_t,
        SA: *mut int32_t,
        n: int32_t,
        k: int32_t,
        fs: int32_t,
    ) -> int32_t;
    fn libsais16_bwt(
        T: *const uint16_t,
        U: *mut uint16_t,
        A: *mut int32_t,
        n: int32_t,
        fs: int32_t,
        freq: *mut int32_t,
    ) -> int32_t;
    fn libsais16_bwt_aux(
        T: *const uint16_t,
        U: *mut uint16_t,
        A: *mut int32_t,
        n: int32_t,
        fs: int32_t,
        freq: *mut int32_t,
        r: int32_t,
        I: *mut int32_t,
    ) -> int32_t;
    fn libsais16_unbwt_aux(
        T: *const uint16_t,
        U: *mut uint16_t,
        A: *mut int32_t,
        n: int32_t,
        freq: *const int32_t,
        r: int32_t,
        I: *const int32_t,
    ) -> int32_t;
    fn malloc(_: std::ffi::c_ulong) -> *mut std::ffi::c_void;
    fn free(_: *mut std::ffi::c_void);
    fn memcpy(
        _: *mut std::ffi::c_void,
        _: *const std::ffi::c_void,
        _: std::ffi::c_ulong,
    ) -> *mut std::ffi::c_void;
    fn memmove(
        _: *mut std::ffi::c_void,
        _: *const std::ffi::c_void,
        _: std::ffi::c_ulong,
    ) -> *mut std::ffi::c_void;
    fn memset(
        _: *mut std::ffi::c_void,
        _: std::ffi::c_int,
        _: std::ffi::c_ulong,
    ) -> *mut std::ffi::c_void;
}
pub type __uint8_t = std::ffi::c_uchar;
pub type __uint16_t = std::ffi::c_ushort;
pub type __int32_t = std::ffi::c_int;
pub type __uint32_t = std::ffi::c_uint;
pub type __int64_t = std::ffi::c_long;
pub type __uint64_t = std::ffi::c_ulong;
pub type int32_t = __int32_t;
pub type int64_t = __int64_t;
pub type uint8_t = __uint8_t;
pub type uint16_t = __uint16_t;
pub type uint32_t = __uint32_t;
pub type uint64_t = __uint64_t;
pub type sa_sint_t = int64_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub union LIBSAIS_THREAD_STATE {
    pub state: C2RustUnnamed,
    pub padding: [uint8_t; 64],
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct C2RustUnnamed {
    pub position: fast_sint_t,
    pub count: fast_sint_t,
    pub m: fast_sint_t,
    pub last_lms_suffix: fast_sint_t,
    pub buckets: *mut sa_sint_t,
    pub cache: *mut LIBSAIS_THREAD_CACHE,
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct LIBSAIS_THREAD_CACHE {
    pub symbol: sa_sint_t,
    pub index: sa_sint_t,
}
pub type fast_sint_t = int64_t;
pub type size_t = std::ffi::c_ulong;
pub type ptrdiff_t = std::ffi::c_long;
pub type sa_uint_t = uint64_t;
pub type fast_uint_t = uint64_t;
unsafe fn libsais16x64_prefetchr(_: *const std::ffi::c_void) {}
unsafe fn libsais16x64_prefetchw(_: *const std::ffi::c_void) {}
unsafe extern "C" fn libsais16x64_align_up(
    mut address: *const std::ffi::c_void,
    mut alignment: size_t,
) -> *mut std::ffi::c_void {
    ((address as ptrdiff_t + alignment as ptrdiff_t
        - 1 as std::ffi::c_int as std::ffi::c_long) & -(alignment as ptrdiff_t))
        as *mut std::ffi::c_void
}
unsafe extern "C" fn libsais16x64_alloc_aligned(
    mut size: size_t,
    mut alignment: size_t,
) -> *mut std::ffi::c_void {
    let mut address: *mut std::ffi::c_void = malloc(
        size
            .wrapping_add(
                ::core::mem::size_of::<std::ffi::c_short>() as std::ffi::c_ulong,
            )
            .wrapping_add(alignment)
            .wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong),
    );
    if !address.is_null() {
        let mut aligned_address: *mut std::ffi::c_void = libsais16x64_align_up(
            (address as ptrdiff_t
                + ::core::mem::size_of::<std::ffi::c_short>() as std::ffi::c_ulong
                    as ptrdiff_t) as *mut std::ffi::c_void,
            alignment,
        );
        *(aligned_address as *mut std::ffi::c_short)
            .offset(
                -(1 as std::ffi::c_int) as isize,
            ) = (aligned_address as ptrdiff_t - address as ptrdiff_t)
            as std::ffi::c_short;
        return aligned_address;
    }
    std::ptr::null_mut::<std::ffi::c_void>()
}
unsafe extern "C" fn libsais16x64_free_aligned(
    mut aligned_address: *mut std::ffi::c_void,
) {
    if !aligned_address.is_null() {
        free(
            (aligned_address as ptrdiff_t
                - *(aligned_address as *mut std::ffi::c_short)
                    .offset(-(1 as std::ffi::c_int) as isize) as std::ffi::c_long)
                as *mut std::ffi::c_void,
        );
    }
}
unsafe extern "C" fn libsais16x64_alloc_thread_state(
    mut threads: sa_sint_t,
) -> *mut LIBSAIS_THREAD_STATE {
    let mut thread_state: *mut LIBSAIS_THREAD_STATE = libsais16x64_alloc_aligned(
        (threads as size_t)
            .wrapping_mul(
                ::core::mem::size_of::<LIBSAIS_THREAD_STATE>() as std::ffi::c_ulong,
            ),
        4096 as std::ffi::c_int as size_t,
    ) as *mut LIBSAIS_THREAD_STATE;
    let mut thread_buckets: *mut sa_sint_t = libsais16x64_alloc_aligned(
        (threads as size_t)
            .wrapping_mul(4 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong,
            )
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        4096 as std::ffi::c_int as size_t,
    ) as *mut sa_sint_t;
    let mut thread_cache: *mut LIBSAIS_THREAD_CACHE = libsais16x64_alloc_aligned(
        (threads as size_t)
            .wrapping_mul(2097184 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(
                ::core::mem::size_of::<LIBSAIS_THREAD_CACHE>() as std::ffi::c_ulong,
            ),
        4096 as std::ffi::c_int as size_t,
    ) as *mut LIBSAIS_THREAD_CACHE;
    if !thread_state.is_null() && !thread_buckets.is_null() && !thread_cache.is_null() {
        let mut t: fast_sint_t = 0;
        t = 0 as std::ffi::c_int as fast_sint_t;
        while t < threads {
            let fresh0 = &mut (*thread_state.offset(t as isize)).state.buckets;
            *fresh0 = thread_buckets;
            thread_buckets = thread_buckets
                .offset(
                    (4 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                );
            let fresh1 = &mut (*thread_state.offset(t as isize)).state.cache;
            *fresh1 = thread_cache;
            thread_cache = thread_cache.offset(2097184 as std::ffi::c_int as isize);
            t += 1;
        }
        return thread_state;
    }
    libsais16x64_free_aligned(thread_cache as *mut std::ffi::c_void);
    libsais16x64_free_aligned(thread_buckets as *mut std::ffi::c_void);
    libsais16x64_free_aligned(thread_state as *mut std::ffi::c_void);
    std::ptr::null_mut::<LIBSAIS_THREAD_STATE>()
}
unsafe extern "C" fn libsais16x64_free_thread_state(
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if !thread_state.is_null() {
        libsais16x64_free_aligned(
            (*thread_state.offset(0 as std::ffi::c_int as isize)).state.cache
                as *mut std::ffi::c_void,
        );
        libsais16x64_free_aligned(
            (*thread_state.offset(0 as std::ffi::c_int as isize)).state.buckets
                as *mut std::ffi::c_void,
        );
        libsais16x64_free_aligned(thread_state as *mut std::ffi::c_void);
    }
}
unsafe extern "C" fn libsais16x64_flip_suffix_markers_omp(
    mut SA: *mut sa_sint_t,
    mut l: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (l / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        l - omp_block_start
    };
    let mut i: fast_sint_t = 0;
    i = omp_block_start;
    while i < omp_block_start + omp_block_size {
        let fresh2 = &mut (*SA.offset(i as isize));
        *fresh2
            ^= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        i += 1;
    }
}
unsafe extern "C" fn libsais16x64_gather_lms_suffixes_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    if omp_block_size > 0 as std::ffi::c_int as std::ffi::c_long {
        let prefetch_distance: fast_sint_t = 128 as std::ffi::c_int as fast_sint_t;
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = omp_block_start + omp_block_size;
        let mut c0: fast_sint_t = *T
            .offset(
                (omp_block_start + omp_block_size
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) as fast_sint_t;
        let mut c1: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
        while j < n
            && {
                c1 = *T.offset(j as isize) as fast_sint_t;
                c1 == c0
            }
        {
            j += 1;
        }
        let mut f0: fast_uint_t = (c0 >= c1) as std::ffi::c_int as fast_uint_t;
        let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
        i = omp_block_start + omp_block_size - 2 as std::ffi::c_int as std::ffi::c_long;
        j = omp_block_start + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchr(
                &*T.offset((i - prefetch_distance) as isize) as *const uint16_t
                    as *const std::ffi::c_void,
            );
            c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            c1 = c0;
            c0 = *T.offset(i as isize) as fast_sint_t;
            f1 = f0;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_gather_lms_suffixes_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_gather_lms_suffixes_16u(
            T,
            SA,
            n,
            n - 1 as std::ffi::c_int as std::ffi::c_long,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_gather_lms_suffixes_32s(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: sa_sint_t = n - 2 as std::ffi::c_int as std::ffi::c_long;
    let mut m: sa_sint_t = n - 1 as std::ffi::c_int as std::ffi::c_long;
    let mut f0: fast_uint_t = 1 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T
        .offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    while i >= 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchr(
            &*T.offset((i - prefetch_distance) as isize) as *const sa_sint_t
                as *const std::ffi::c_void,
        );
        c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m -= (f1 & !f0) as sa_sint_t;
        c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
        m -= (f0 & !f1) as sa_sint_t;
        c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
        m -= (f1 & !f0) as sa_sint_t;
        c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
        m -= (f0 & !f1) as sa_sint_t;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        c1 = c0;
        c0 = *T.offset(i as isize);
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m -= (f0 & !f1) as sa_sint_t;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    n - 1 as std::ffi::c_int as std::ffi::c_long - m
}
unsafe extern "C" fn libsais16x64_gather_compacted_lms_suffixes_32s(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: sa_sint_t = n - 2 as std::ffi::c_int as std::ffi::c_long;
    let mut m: sa_sint_t = n - 1 as std::ffi::c_int as std::ffi::c_long;
    let mut f0: fast_uint_t = 1 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T
        .offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    while i >= 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchr(
            &*T.offset((i - prefetch_distance) as isize) as *const sa_sint_t
                as *const std::ffi::c_void,
        );
        c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m
            -= (f1 & !f0
                & (c0 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                    as std::ffi::c_ulong) as sa_sint_t;
        c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
        m
            -= (f0 & !f1
                & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                    as std::ffi::c_ulong) as sa_sint_t;
        c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
        m
            -= (f1 & !f0
                & (c0 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                    as std::ffi::c_ulong) as sa_sint_t;
        c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
        m
            -= (f0 & !f1
                & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                    as std::ffi::c_ulong) as sa_sint_t;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        c1 = c0;
        c0 = *T.offset(i as isize);
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m
            -= (f0 & !f1
                & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                    as std::ffi::c_ulong) as sa_sint_t;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    n - 1 as std::ffi::c_int as std::ffi::c_long - m
}
unsafe extern "C" fn libsais16x64_count_lms_suffixes_32s_2k(
    mut T: *const sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (2 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut i: sa_sint_t = n - 2 as std::ffi::c_int as std::ffi::c_long;
    let mut f0: fast_uint_t = 1 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T
        .offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    while i >= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchr(
            &*T
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    ((*T
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    ((*T
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    ((*T
                        .offset(
                            (i - prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    ((*T
                        .offset(
                            (i - prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh3 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f1 & !f0) as fast_sint_t) as isize,
            ));
        *fresh3 += 1;
        c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh4 = &mut (*buckets
            .offset(
                (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f0 & !f1) as fast_sint_t) as isize,
            ));
        *fresh4 += 1;
        c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh5 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f1 & !f0) as fast_sint_t) as isize,
            ));
        *fresh5 += 1;
        c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh6 = &mut (*buckets
            .offset(
                (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f0 & !f1) as fast_sint_t) as isize,
            ));
        *fresh6 += 1;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        c1 = c0;
        c0 = *T.offset(i as isize);
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh7 = &mut (*buckets
            .offset(
                (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f0 & !f1) as fast_sint_t) as isize,
            ));
        *fresh7 += 1;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    let fresh8 = &mut (*buckets
        .offset(
            (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh8 += 1;
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (4 as std::ffi::c_int as size_t)
            .wrapping_mul(
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong,
            )
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut m: fast_sint_t = omp_block_start + omp_block_size
        - 1 as std::ffi::c_int as std::ffi::c_long;
    if omp_block_size > 0 as std::ffi::c_int as std::ffi::c_long {
        let prefetch_distance: fast_sint_t = 128 as std::ffi::c_int as fast_sint_t;
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = m + 1 as std::ffi::c_int as std::ffi::c_long;
        let mut c0: fast_sint_t = *T.offset(m as isize) as fast_sint_t;
        let mut c1: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
        while j < n
            && {
                c1 = *T.offset(j as isize) as fast_sint_t;
                c1 == c0
            }
        {
            j += 1;
        }
        let mut f0: fast_uint_t = (c0 >= c1) as std::ffi::c_int as fast_uint_t;
        let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
        i = m - 1 as std::ffi::c_int as std::ffi::c_long;
        j = omp_block_start + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchr(
                &*T.offset((i - prefetch_distance) as isize) as *const uint16_t
                    as *const std::ffi::c_void,
            );
            c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh9 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
                ));
            *fresh9 += 1;
            c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh10 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh10 += 1;
            c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh11 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
                ));
            *fresh11 += 1;
            c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh12 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh12 += 1;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            c1 = c0;
            c0 = *T.offset(i as isize) as fast_sint_t;
            f1 = f0;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh13 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh13 += 1;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c1 = (if i >= 0 as std::ffi::c_int as std::ffi::c_long {
            *T.offset(i as isize) as std::ffi::c_int
        } else {
            -(1 as std::ffi::c_int)
        }) as fast_sint_t;
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
            as fast_sint_t;
        let fresh14 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                    + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
            ));
        *fresh14 += 1;
    }
    omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long
        - m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        m = libsais16x64_count_and_gather_lms_suffixes_16u(
            T,
            SA,
            n,
            buckets,
            omp_block_start,
            omp_block_size,
        );
    }
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_4k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (4 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut m: fast_sint_t = omp_block_start + omp_block_size
        - 1 as std::ffi::c_int as std::ffi::c_long;
    if omp_block_size > 0 as std::ffi::c_int as std::ffi::c_long {
        let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = m + 1 as std::ffi::c_int as std::ffi::c_long;
        let mut c0: fast_sint_t = *T.offset(m as isize);
        let mut c1: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
        while j < n
            && {
                c1 = *T.offset(j as isize);
                c1 == c0
            }
        {
            j += 1;
        }
        let mut f0: fast_uint_t = (c0 >= c1) as std::ffi::c_int as fast_uint_t;
        let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
        i = m - 1 as std::ffi::c_int as std::ffi::c_long;
        j = omp_block_start + prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchr(
                &*T
                    .offset(
                        (i
                            - 2 as std::ffi::c_int as std::ffi::c_long
                                * prefetch_distance) as isize,
                    ) as *const sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh15 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
                ));
            *fresh15 += 1;
            c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh16 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh16 += 1;
            c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh17 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
                ));
            *fresh17 += 1;
            c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh18 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh18 += 1;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            c1 = c0;
            c0 = *T.offset(i as isize);
            f1 = f0;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh19 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                        + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
                ));
            *fresh19 += 1;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c1 = if i >= 0 as std::ffi::c_int as std::ffi::c_long {
            *T.offset(i as isize)
        } else {
            -(1 as std::ffi::c_int) as std::ffi::c_long
        };
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
            as fast_sint_t;
        let fresh20 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                    + f0.wrapping_add(f0).wrapping_add(f1) as fast_sint_t) as isize,
            ));
        *fresh20 += 1;
    }
    omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long
        - m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_2k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (2 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut m: fast_sint_t = omp_block_start + omp_block_size
        - 1 as std::ffi::c_int as std::ffi::c_long;
    if omp_block_size > 0 as std::ffi::c_int as std::ffi::c_long {
        let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = m + 1 as std::ffi::c_int as std::ffi::c_long;
        let mut c0: fast_sint_t = *T.offset(m as isize);
        let mut c1: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
        while j < n
            && {
                c1 = *T.offset(j as isize);
                c1 == c0
            }
        {
            j += 1;
        }
        let mut f0: fast_uint_t = (c0 >= c1) as std::ffi::c_int as fast_uint_t;
        let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
        i = m - 1 as std::ffi::c_int as std::ffi::c_long;
        j = omp_block_start + prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchr(
                &*T
                    .offset(
                        (i
                            - 2 as std::ffi::c_int as std::ffi::c_long
                                * prefetch_distance) as isize,
                    ) as *const sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        ((*T
                            .offset(
                                (i - prefetch_distance
                                    - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh21 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f1 & !f0) as fast_sint_t) as isize,
                ));
            *fresh21 += 1;
            c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh22 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh22 += 1;
            c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
                as fast_sint_t;
            let fresh23 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f1 & !f0) as fast_sint_t) as isize,
                ));
            *fresh23 += 1;
            c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh24 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh24 += 1;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            c1 = c0;
            c0 = *T.offset(i as isize);
            f1 = f0;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong).wrapping_sub(f0 & !f1) as fast_sint_t
                as fast_sint_t;
            let fresh25 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh25 += 1;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c1 = if i >= 0 as std::ffi::c_int as std::ffi::c_long {
            *T.offset(i as isize)
        } else {
            -(1 as std::ffi::c_int) as std::ffi::c_long
        };
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m = (m as std::ffi::c_ulong).wrapping_sub(f1 & !f0) as fast_sint_t
            as fast_sint_t;
        let fresh26 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f1 & !f0) as fast_sint_t) as isize,
            ));
        *fresh26 += 1;
    }
    omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long
        - m
}
unsafe extern "C" fn libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (2 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut m: fast_sint_t = omp_block_start + omp_block_size
        - 1 as std::ffi::c_int as std::ffi::c_long;
    if omp_block_size > 0 as std::ffi::c_int as std::ffi::c_long {
        let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = m + 1 as std::ffi::c_int as std::ffi::c_long;
        let mut c0: fast_sint_t = *T.offset(m as isize);
        let mut c1: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
        while j < n
            && {
                c1 = *T.offset(j as isize);
                c1 == c0
            }
        {
            j += 1;
        }
        let mut f0: fast_uint_t = (c0 >= c1) as std::ffi::c_int as fast_uint_t;
        let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
        i = m - 1 as std::ffi::c_int as std::ffi::c_long;
        j = omp_block_start + prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchr(
                &*T
                    .offset(
                        (i
                            - 2 as std::ffi::c_int as std::ffi::c_long
                                * prefetch_distance) as isize,
                    ) as *const sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        (((*T
                            .offset(
                                (i - prefetch_distance
                                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) & 9223372036854775807 as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        (((*T
                            .offset(
                                (i - prefetch_distance
                                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) & 9223372036854775807 as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        (((*T
                            .offset(
                                (i - prefetch_distance
                                    - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) & 9223372036854775807 as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *buckets
                    .offset(
                        (((*T
                            .offset(
                                (i - prefetch_distance
                                    - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) & 9223372036854775807 as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong)
                .wrapping_sub(
                    f1 & !f0
                        & (c0 >= 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_ulong,
                ) as fast_sint_t as fast_sint_t;
            c0 &= 9223372036854775807 as std::ffi::c_long;
            let fresh27 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f1 & !f0) as fast_sint_t) as isize,
                ));
            *fresh27 += 1;
            c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong)
                .wrapping_sub(
                    f0 & !f1
                        & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_ulong,
                ) as fast_sint_t as fast_sint_t;
            c1 &= 9223372036854775807 as std::ffi::c_long;
            let fresh28 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh28 += 1;
            c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
            f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong)
                .wrapping_sub(
                    f1 & !f0
                        & (c0 >= 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_ulong,
                ) as fast_sint_t as fast_sint_t;
            c0 &= 9223372036854775807 as std::ffi::c_long;
            let fresh29 = &mut (*buckets
                .offset(
                    (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f1 & !f0) as fast_sint_t) as isize,
                ));
            *fresh29 += 1;
            c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong)
                .wrapping_sub(
                    f0 & !f1
                        & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_ulong,
                ) as fast_sint_t as fast_sint_t;
            c1 &= 9223372036854775807 as std::ffi::c_long;
            let fresh30 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh30 += 1;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            c1 = c0;
            c0 = *T.offset(i as isize);
            f1 = f0;
            f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
            *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m = (m as std::ffi::c_ulong)
                .wrapping_sub(
                    f0 & !f1
                        & (c1 >= 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_ulong,
                ) as fast_sint_t as fast_sint_t;
            c1 &= 9223372036854775807 as std::ffi::c_long;
            let fresh31 = &mut (*buckets
                .offset(
                    (((c1 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                        + (f0 & !f1) as fast_sint_t) as isize,
                ));
            *fresh31 += 1;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c1 = if i >= 0 as std::ffi::c_int as std::ffi::c_long {
            *T.offset(i as isize)
        } else {
            -(1 as std::ffi::c_int) as std::ffi::c_long
        };
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        *SA.offset(m as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
        m = (m as std::ffi::c_ulong)
            .wrapping_sub(
                f1 & !f0
                    & (c0 >= 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                        as std::ffi::c_ulong,
            ) as fast_sint_t as fast_sint_t;
        c0 &= 9223372036854775807 as std::ffi::c_long;
        let fresh32 = &mut (*buckets
            .offset(
                (((c0 as fast_uint_t as fast_sint_t) << 1 as std::ffi::c_int)
                    + (f1 & !f0) as fast_sint_t) as isize,
            ));
        *fresh32 += 1;
    }
    omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long
        - m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_4k_nofs_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        m = libsais16x64_count_and_gather_lms_suffixes_32s_4k(
            T,
            SA,
            n,
            k,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_2k_nofs_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        m = libsais16x64_count_and_gather_lms_suffixes_32s_2k(
            T,
            SA,
            n,
            k,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        m = libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k(
            T,
            SA,
            n,
            k,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_4k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0;
    m = libsais16x64_count_and_gather_lms_suffixes_32s_4k_nofs_omp(
        T,
        SA,
        n,
        k,
        buckets,
        threads,
    );
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_lms_suffixes_32s_2k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut m: sa_sint_t = 0;
    m = libsais16x64_count_and_gather_lms_suffixes_32s_2k_nofs_omp(
        T,
        SA,
        n,
        k,
        buckets,
        threads,
    );
    m
}
unsafe extern "C" fn libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_nofs_omp(
        T,
        SA,
        n,
        k,
        buckets,
        threads,
    );
}
unsafe extern "C" fn libsais16x64_count_suffixes_32s(
    mut T: *const sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = 0 as std::ffi::c_int as fast_sint_t;
    j = n - 7 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &*T.offset((i + prefetch_distance) as isize) as *const sa_sint_t
                as *const std::ffi::c_void,
        );
        let fresh33 = &mut (*buckets
            .offset(
                *T.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh33 += 1;
        let fresh34 = &mut (*buckets
            .offset(
                *T.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh34 += 1;
        let fresh35 = &mut (*buckets
            .offset(
                *T.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh35 += 1;
        let fresh36 = &mut (*buckets
            .offset(
                *T.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh36 += 1;
        let fresh37 = &mut (*buckets
            .offset(
                *T.offset((i + 4 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh37 += 1;
        let fresh38 = &mut (*buckets
            .offset(
                *T.offset((i + 5 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh38 += 1;
        let fresh39 = &mut (*buckets
            .offset(
                *T.offset((i + 6 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh39 += 1;
        let fresh40 = &mut (*buckets
            .offset(
                *T.offset((i + 7 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh40 += 1;
        i += 8 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 7 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let fresh41 = &mut (*buckets.offset(*T.offset(i as isize) as isize));
        *fresh41 += 1;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_start_and_end_16u(
    mut buckets: *mut sa_sint_t,
    mut freq: *mut sa_sint_t,
) -> sa_sint_t {
    let mut bucket_start: *mut sa_sint_t = &mut *buckets
        .offset(
            (6 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut bucket_end: *mut sa_sint_t = &mut *buckets
        .offset(
            (7 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut k: fast_sint_t = -(1 as std::ffi::c_int) as fast_sint_t;
    if !freq.is_null() {
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = 0;
        let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
        i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t;
        j = 0 as std::ffi::c_int as fast_sint_t;
        while i
            <= (((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                << 8 as std::ffi::c_int) as fast_sint_t
                - 1 as std::ffi::c_int as std::ffi::c_long) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t
        {
            let mut total: sa_sint_t = *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 2 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 3 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
            *bucket_start.offset(j as isize) = sum;
            sum += total;
            *bucket_end.offset(j as isize) = sum;
            k = if total > 0 as std::ffi::c_int as std::ffi::c_long { j } else { k };
            *freq.offset(j as isize) = total;
            i
                += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t;
            j += 1 as std::ffi::c_int as std::ffi::c_long;
        }
    } else {
        let mut i_0: fast_sint_t = 0;
        let mut j_0: fast_sint_t = 0;
        let mut sum_0: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
        i_0 = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t;
        j_0 = 0 as std::ffi::c_int as fast_sint_t;
        while i_0
            <= (((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                << 8 as std::ffi::c_int) as fast_sint_t
                - 1 as std::ffi::c_int as std::ffi::c_long) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t
        {
            let mut total_0: sa_sint_t = *buckets
                .offset(
                    (i_0
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i_0
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i_0
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 2 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i_0
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 3 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
            *bucket_start.offset(j_0 as isize) = sum_0;
            sum_0 += total_0;
            *bucket_end.offset(j_0 as isize) = sum_0;
            k = if total_0 > 0 as std::ffi::c_int as std::ffi::c_long { j_0 } else { k };
            i_0
                += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t;
            j_0 += 1 as std::ffi::c_int as std::ffi::c_long;
        }
    }
    k + 1 as std::ffi::c_int as std::ffi::c_long
}
unsafe extern "C" fn libsais16x64_initialize_buckets_start_and_end_32s_6k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut bucket_start: *mut sa_sint_t = &mut *buckets
        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut bucket_end: *mut sa_sint_t = &mut *buckets
        .offset((5 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        *bucket_start.offset(j as isize) = sum;
        sum
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 2 as std::ffi::c_int as fast_sint_t)) as isize,
                    )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 3 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *bucket_end.offset(j as isize) = sum;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_start_and_end_32s_4k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut bucket_start: *mut sa_sint_t = &mut *buckets
        .offset((2 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut bucket_end: *mut sa_sint_t = &mut *buckets
        .offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        *bucket_start.offset(j as isize) = sum;
        sum
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 1 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *bucket_end.offset(j as isize) = sum;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_end_32s_2k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut sum0: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        sum0
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 1 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_start_and_end_32s_2k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        *buckets.offset(j as isize) = *buckets.offset(i as isize);
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    *buckets.offset(k as isize) = 0 as std::ffi::c_int as sa_sint_t;
    memcpy(
        &mut *buckets.offset((k + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as *mut sa_sint_t as *mut std::ffi::c_void,
        buckets as *const std::ffi::c_void,
        (k as size_t)
            .wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_initialize_buckets_start_32s_1k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = 0 as std::ffi::c_int as fast_sint_t;
    while i <= k - 1 as std::ffi::c_int as std::ffi::c_long {
        let mut tmp: sa_sint_t = *buckets.offset(i as isize);
        *buckets.offset(i as isize) = sum;
        sum += tmp;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_end_32s_1k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = 0 as std::ffi::c_int as fast_sint_t;
    while i <= k - 1 as std::ffi::c_int as std::ffi::c_long {
        sum += *buckets.offset(i as isize);
        *buckets.offset(i as isize) = sum;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_16u(
    mut T: *const uint16_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
) -> sa_sint_t {
    let mut f0: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T.offset(first_lms_suffix as isize) as fast_sint_t;
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    loop {
        first_lms_suffix -= 1;
        if first_lms_suffix < 0 as std::ffi::c_int as std::ffi::c_long {
            break;
        }
        c1 = c0;
        c0 = *T.offset(first_lms_suffix as isize) as fast_sint_t;
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh42 = &mut (*buckets
            .offset(
                (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                    + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
            ));
        *fresh42 -= 1;
    }
    let fresh43 = &mut (*buckets
        .offset(
            (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                + f0.wrapping_add(f0) as fast_sint_t) as isize,
        ));
    *fresh43 -= 1;
    let mut temp_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (4 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= (((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
            as fast_sint_t - 1 as std::ffi::c_int as std::ffi::c_long)
            << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t
    {
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum;
        sum
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 3 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *temp_bucket.offset(j as isize) = sum;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
    sum
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(
    mut T: *const sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
) {
    let fresh44 = &mut (*buckets
        .offset(
            ((*T.offset(first_lms_suffix as isize) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh44 += 1;
    let fresh45 = &mut (*buckets
        .offset(
            ((*T.offset(first_lms_suffix as isize) << 1 as std::ffi::c_int)
                + 1 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh45 -= 1;
    let mut i: fast_sint_t = 0;
    let mut sum0: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut sum1: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        sum0
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 1 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        sum1
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                );
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum1;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(
    mut T: *const sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
) -> sa_sint_t {
    let mut f0: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T.offset(first_lms_suffix as isize);
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    loop {
        first_lms_suffix -= 1;
        if first_lms_suffix < 0 as std::ffi::c_int as std::ffi::c_long {
            break;
        }
        c1 = c0;
        c0 = *T.offset(first_lms_suffix as isize);
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        let fresh46 = &mut (*buckets
            .offset(
                (((c1 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                    + f1.wrapping_add(f1).wrapping_add(f0) as fast_sint_t) as isize,
            ));
        *fresh46 -= 1;
    }
    let fresh47 = &mut (*buckets
        .offset(
            (((c0 as fast_uint_t as fast_sint_t) << 2 as std::ffi::c_int)
                + f0.wrapping_add(f0) as fast_sint_t) as isize,
        ));
    *fresh47 -= 1;
    let mut temp_bucket: *mut sa_sint_t = &mut *buckets
        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        sum
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 3 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *temp_bucket.offset(j as isize) = sum;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    sum
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_radix_and_partial_sorting_32s_4k(
    mut T: *const sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
) {
    let mut bucket_start: *mut sa_sint_t = &mut *buckets
        .offset((2 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut bucket_end: *mut sa_sint_t = &mut *buckets
        .offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let fresh48 = &mut (*buckets
        .offset(
            ((*T.offset(first_lms_suffix as isize) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh48 += 1;
    let fresh49 = &mut (*buckets
        .offset(
            ((*T.offset(first_lms_suffix as isize) << 1 as std::ffi::c_int)
                + 1 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh49 -= 1;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum0: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut sum1: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        *bucket_start.offset(j as isize) = sum1;
        sum0
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                );
        sum1
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 1 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *bucket_end.offset(j as isize) = sum1;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const uint16_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const uint16_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const uint16_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const uint16_t as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh50 = &mut (*induction_bucket
            .offset(
                (((*T.offset(p0 as isize) as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh50 -= 1;
        *SA.offset(*fresh50 as isize) = p0;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh51 = &mut (*induction_bucket
            .offset(
                (((*T.offset(p1 as isize) as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh51 -= 1;
        *SA.offset(*fresh51 as isize) = p1;
        let mut p2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh52 = &mut (*induction_bucket
            .offset(
                (((*T.offset(p2 as isize) as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh52 -= 1;
        *SA.offset(*fresh52 as isize) = p2;
        let mut p3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh53 = &mut (*induction_bucket
            .offset(
                (((*T.offset(p3 as isize) as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh53 -= 1;
        *SA.offset(*fresh53 as isize) = p3;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        let fresh54 = &mut (*induction_bucket
            .offset(
                (((*T.offset(p as isize) as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh54 -= 1;
        *SA.offset(*fresh54 as isize) = p;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut flags: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
        let fresh55 = &mut (*buckets
            .offset(
                (4 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as isize,
            ));
        *fresh55 -= 1;
    }
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_radix_sort_lms_suffixes_16u(
            T,
            SA,
            &mut *buckets
                .offset(
                    (4 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            n - m + 1 as std::ffi::c_int as std::ffi::c_long,
            m - 1 as std::ffi::c_int as std::ffi::c_long,
        );
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_32s_6k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    *T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    *T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    *T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    *T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh56 = &mut (*induction_bucket.offset(*T.offset(p0 as isize) as isize));
        *fresh56 -= 1;
        *SA.offset(*fresh56 as isize) = p0;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh57 = &mut (*induction_bucket.offset(*T.offset(p1 as isize) as isize));
        *fresh57 -= 1;
        *SA.offset(*fresh57 as isize) = p1;
        let mut p2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh58 = &mut (*induction_bucket.offset(*T.offset(p2 as isize) as isize));
        *fresh58 -= 1;
        *SA.offset(*fresh58 as isize) = p2;
        let mut p3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh59 = &mut (*induction_bucket.offset(*T.offset(p3 as isize) as isize));
        *fresh59 -= 1;
        *SA.offset(*fresh59 as isize) = p3;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        let fresh60 = &mut (*induction_bucket.offset(*T.offset(p as isize) as isize));
        *fresh60 -= 1;
        *SA.offset(*fresh60 as isize) = p;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_32s_2k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    ((*T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    ((*T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    ((*T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *induction_bucket
                .offset(
                    ((*T
                        .offset(
                            *SA
                                .offset(
                                    (i - prefetch_distance
                                        - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                                ) as isize,
                        ) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh61 = &mut (*induction_bucket
            .offset(
                ((*T.offset(p0 as isize) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh61 -= 1;
        *SA.offset(*fresh61 as isize) = p0;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh62 = &mut (*induction_bucket
            .offset(
                ((*T.offset(p1 as isize) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh62 -= 1;
        *SA.offset(*fresh62 as isize) = p1;
        let mut p2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh63 = &mut (*induction_bucket
            .offset(
                ((*T.offset(p2 as isize) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh63 -= 1;
        *SA.offset(*fresh63 as isize) = p2;
        let mut p3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh64 = &mut (*induction_bucket
            .offset(
                ((*T.offset(p3 as isize) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh64 -= 1;
        *SA.offset(*fresh64 as isize) = p3;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        let fresh65 = &mut (*induction_bucket
            .offset(
                ((*T.offset(p as isize) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ));
        *fresh65 -= 1;
        *SA.offset(*fresh65 as isize) = p;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_32s_6k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || m < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_radix_sort_lms_suffixes_32s_6k(
            T,
            SA,
            induction_bucket,
            n - m + 1 as std::ffi::c_int as std::ffi::c_long,
            m - 1 as std::ffi::c_int as std::ffi::c_long,
        );
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_32s_2k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || m < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_radix_sort_lms_suffixes_32s_2k(
            T,
            SA,
            induction_bucket,
            n - m + 1 as std::ffi::c_int as std::ffi::c_long,
            m - 1 as std::ffi::c_int as std::ffi::c_long,
        );
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_lms_suffixes_32s_1k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: sa_sint_t = n - 2 as std::ffi::c_int as std::ffi::c_long;
    let mut m: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut f0: fast_uint_t = 1 as std::ffi::c_int as fast_uint_t;
    let mut f1: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    let mut c0: fast_sint_t = *T
        .offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
    let mut c1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut c2: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    while i >= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchr(
            &*T
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    *T
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    *T
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    *T
                        .offset(
                            (i - prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *buckets
                .offset(
                    *T
                        .offset(
                            (i - prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        c1 = *T.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        if f1 & !f0 != 0 {
            c2 = c0;
            let fresh66 = &mut (*buckets.offset(c2 as isize));
            *fresh66 -= 1;
            *SA.offset(*fresh66 as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m += 1;
        }
        c0 = *T.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        if f0 & !f1 != 0 {
            c2 = c1;
            let fresh67 = &mut (*buckets.offset(c2 as isize));
            *fresh67 -= 1;
            *SA.offset(*fresh67 as isize) = i - 0 as std::ffi::c_int as std::ffi::c_long;
            m += 1;
        }
        c1 = *T.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        f1 = (c1 > c0 - f0 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        if f1 & !f0 != 0 {
            c2 = c0;
            let fresh68 = &mut (*buckets.offset(c2 as isize));
            *fresh68 -= 1;
            *SA.offset(*fresh68 as isize) = i - 1 as std::ffi::c_int as std::ffi::c_long;
            m += 1;
        }
        c0 = *T.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        if f0 & !f1 != 0 {
            c2 = c1;
            let fresh69 = &mut (*buckets.offset(c2 as isize));
            *fresh69 -= 1;
            *SA.offset(*fresh69 as isize) = i - 2 as std::ffi::c_int as std::ffi::c_long;
            m += 1;
        }
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        c1 = c0;
        c0 = *T.offset(i as isize);
        f1 = f0;
        f0 = (c0 > c1 - f1 as fast_sint_t) as std::ffi::c_int as fast_uint_t;
        if f0 & !f1 != 0 {
            c2 = c1;
            let fresh70 = &mut (*buckets.offset(c2 as isize));
            *fresh70 -= 1;
            *SA.offset(*fresh70 as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            m += 1;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    if m > 1 as std::ffi::c_int as std::ffi::c_long {
        *SA
            .offset(
                *buckets.offset(c2 as isize) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
    }
    m
}
unsafe extern "C" fn libsais16x64_radix_sort_set_markers_32s_6k(
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *induction_bucket
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let fresh71 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh71
            |= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        let fresh72 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh72
            |= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        let fresh73 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh73
            |= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        let fresh74 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ));
        *fresh74
            |= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let fresh75 = &mut (*SA.offset(*induction_bucket.offset(i as isize) as isize));
        *fresh75
            |= -(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_set_markers_32s_4k(
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *induction_bucket
                .offset(
                    (((i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance) << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (((i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long)
                                << 1 as std::ffi::c_int)
                                + 0 as std::ffi::c_int as fast_sint_t) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (((i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long)
                                << 1 as std::ffi::c_int)
                                + 0 as std::ffi::c_int as fast_sint_t) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (((i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long)
                                << 1 as std::ffi::c_int)
                                + 0 as std::ffi::c_int as fast_sint_t) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    *induction_bucket
                        .offset(
                            (((i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long)
                                << 1 as std::ffi::c_int)
                                + 0 as std::ffi::c_int as fast_sint_t) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let fresh76 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset(
                        (((i + 0 as std::ffi::c_int as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as isize,
            ));
        *fresh76
            |= (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        let fresh77 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset(
                        (((i + 1 as std::ffi::c_int as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as isize,
            ));
        *fresh77
            |= (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        let fresh78 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset(
                        (((i + 2 as std::ffi::c_int as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as isize,
            ));
        *fresh78
            |= (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        let fresh79 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset(
                        (((i + 3 as std::ffi::c_int as std::ffi::c_long)
                            << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as isize,
            ));
        *fresh79
            |= (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let fresh80 = &mut (*SA
            .offset(
                *induction_bucket
                    .offset(
                        ((i << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as isize,
            ));
        *fresh80
            |= (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_radix_sort_set_markers_32s_6k_omp(
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = k - 1 as std::ffi::c_int as std::ffi::c_long;
    libsais16x64_radix_sort_set_markers_32s_6k(
        SA,
        induction_bucket,
        omp_block_start,
        omp_block_size,
    );
}
unsafe extern "C" fn libsais16x64_radix_sort_set_markers_32s_4k_omp(
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = k - 1 as std::ffi::c_int as std::ffi::c_long;
    libsais16x64_radix_sort_set_markers_32s_4k(
        SA,
        induction_bucket,
        omp_block_start,
        omp_block_size,
    );
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_partial_sorting_16u(
    mut T: *const uint16_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
) {
    let mut temp_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (4 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let fresh81 = &mut (*buckets
        .offset(
            (((*T.offset(first_lms_suffix as isize) as fast_uint_t as fast_sint_t)
                << 2 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t) as isize,
        ));
    *fresh81 += 1;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum0: sa_sint_t = left_suffixes_count
        + 1 as std::ffi::c_int as std::ffi::c_long;
    let mut sum1: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= (((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
            as fast_sint_t - 1 as std::ffi::c_int as std::ffi::c_long)
            << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t
    {
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        sum0
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                + *buckets
                    .offset(
                        (i
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 2 as std::ffi::c_int)
                                + 2 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
        sum1
            += *buckets
                .offset(
                    (i
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 2 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                );
        *buckets
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *buckets
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum1;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_initialize_buckets_for_partial_sorting_32s_6k(
    mut T: *const sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
) {
    let mut temp_bucket: *mut sa_sint_t = &mut *buckets
        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut sum0: sa_sint_t = left_suffixes_count
        + 1 as std::ffi::c_int as std::ffi::c_long;
    let mut sum1: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut sum2: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    first_lms_suffix = *T.offset(first_lms_suffix as isize);
    i = ((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    j = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        != (first_lms_suffix << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        let mut SS: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut LS: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut SL: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 2 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut LL: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 3 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum2;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 2 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 3 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
        sum0 += SS + SL;
        sum1 += LS;
        sum2 += LS + LL;
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum1;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
    sum1 += 1 as std::ffi::c_int as std::ffi::c_long;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 2 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        let mut SS_0: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut LS_0: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut SL_0: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 2 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        let mut LL_0: sa_sint_t = *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 3 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum2;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 2 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
        *buckets
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 3 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
        sum0 += SS_0 + SL_0;
        sum1 += LS_0;
        sum2 += LS_0 + LL_0;
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum0;
        *temp_bucket
            .offset(
                (j
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = sum1;
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
        j
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (4 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset(
            (2 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p0 &= 9223372036854775807 as std::ffi::c_long;
        let mut v0: sa_sint_t = ((*T
            .offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                >= *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh82 = &mut (*induction_bucket.offset(v0 as isize));
        let fresh83 = *fresh82;
        *fresh82 += 1;
        *SA
            .offset(
                fresh83 as isize,
            ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v0 as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v0 as isize) = d;
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p1 &= 9223372036854775807 as std::ffi::c_long;
        let mut v1: sa_sint_t = ((*T
            .offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                >= *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh84 = &mut (*induction_bucket.offset(v1 as isize));
        let fresh85 = *fresh84;
        *fresh84 += 1;
        *SA
            .offset(
                fresh85 as isize,
            ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v1 as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v1 as isize) = d;
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        d
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p &= 9223372036854775807 as std::ffi::c_long;
        let mut v: sa_sint_t = ((*T
            .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                >= *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh86 = &mut (*induction_bucket.offset(v as isize));
        let fresh87 = *fresh86;
        *fresh86 += 1;
        *SA
            .offset(
                fresh87 as isize,
            ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v as isize) = d;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (4 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset(
            (2 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let fresh88 = &mut (*induction_bucket
        .offset(
            (((*T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t) << 1 as std::ffi::c_int)
                + (*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int
                    >= *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as std::ffi::c_int) as std::ffi::c_int as fast_sint_t) as isize,
        ));
    let fresh89 = *fresh88;
    *fresh88 += 1;
    *SA
        .offset(
            fresh89 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | (-(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long);
    d += 1;
    *distinct_names
        .offset(
            (((*T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as fast_sint_t) << 1 as std::ffi::c_int)
                + (*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int
                    >= *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as std::ffi::c_int) as std::ffi::c_int as fast_sint_t) as isize,
        ) = d;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || left_suffixes_count < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        d = libsais16x64_partial_sorting_scan_left_to_right_16u(
            T,
            SA,
            buckets,
            d,
            0 as std::ffi::c_int as fast_sint_t,
            left_suffixes_count,
        );
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_6k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size
        - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) & 9223372036854775807 as std::ffi::c_long;
        let mut v0: sa_sint_t = (*T
            .offset(
                (p0
                    - (p0 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                        as std::ffi::c_long) as isize,
            ) << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t;
        libsais16x64_prefetchw(
            &mut *buckets.offset(v0 as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p1: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) & 9223372036854775807 as std::ffi::c_long;
        let mut v1: sa_sint_t = (*T
            .offset(
                (p1
                    - (p1 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                        as std::ffi::c_long) as isize,
            ) << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t;
        libsais16x64_prefetchw(
            &mut *buckets.offset(v1 as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p2: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p2 &= 9223372036854775807 as std::ffi::c_long;
        let mut v2: sa_sint_t = (*T
            .offset((p2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p2 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                >= *T.offset((p2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh90 = &mut (*buckets.offset(v2 as isize));
        let fresh91 = *fresh90;
        *fresh90 += 1;
        *SA
            .offset(
                fresh91 as isize,
            ) = (p2 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets
                .offset((2 as std::ffi::c_int as std::ffi::c_long + v2) as isize) != d)
                as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v2) as isize) = d;
        let mut p3: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p3 &= 9223372036854775807 as std::ffi::c_long;
        let mut v3: sa_sint_t = (*T
            .offset((p3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p3 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                >= *T.offset((p3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh92 = &mut (*buckets.offset(v3 as isize));
        let fresh93 = *fresh92;
        *fresh92 += 1;
        *SA
            .offset(
                fresh93 as isize,
            ) = (p3 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets
                .offset((2 as std::ffi::c_int as std::ffi::c_long + v3) as isize) != d)
                as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v3) as isize) = d;
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        += 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        d
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p &= 9223372036854775807 as std::ffi::c_long;
        let mut v: sa_sint_t = (*T
            .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                >= *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh94 = &mut (*buckets.offset(v as isize));
        let fresh95 = *fresh94;
        *fresh94 += 1;
        *SA
            .offset(
                fresh95 as isize,
            ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v) as isize)
                != d) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v) as isize) = d;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_4k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset((2 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset((0 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size
        - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int))
                } else {
                    2 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            Ts0.offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int))
                } else {
                    2 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            Ts1.offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            let Ts2: fast_sint_t = *T
                .offset(
                    ((s2
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int)))
                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            libsais16x64_prefetchw(
                &mut *induction_bucket.offset(Ts2 as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *distinct_names
                    .offset(
                        ((Ts2 << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            let Ts3: fast_sint_t = *T
                .offset(
                    ((s3
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int)))
                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            libsais16x64_prefetchw(
                &mut *induction_bucket.offset(Ts3 as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *distinct_names
                    .offset(
                        ((Ts3 << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p0 >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p0
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v0: sa_sint_t = (*T
                .offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as fast_sint_t;
            let fresh96 = &mut (*induction_bucket
                .offset(
                    *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh97 = *fresh96;
            *fresh96 += 1;
            *SA
                .offset(
                    fresh97 as isize,
                ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v0 as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v0 as isize) = d;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p1 >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p1
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v1: sa_sint_t = (*T
                .offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as fast_sint_t;
            let fresh98 = &mut (*induction_bucket
                .offset(
                    *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh99 = *fresh98;
            *fresh98 += 1;
            *SA
                .offset(
                    fresh99 as isize,
                ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v1 as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v1 as isize) = d;
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        += 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(i as isize) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v: sa_sint_t = (*T
                .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as fast_sint_t;
            let fresh100 = &mut (*induction_bucket
                .offset(
                    *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh101 = *fresh100;
            *fresh100 += 1;
            *SA
                .offset(
                    fresh101 as isize,
                ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v as isize) = d;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_1k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size
        - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s2 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s3 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh102 = &mut (*induction_bucket
                .offset(
                    *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh103 = *fresh102;
            *fresh102 += 1;
            *SA
                .offset(
                    fresh103 as isize,
                ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh104 = &mut (*induction_bucket
                .offset(
                    *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh105 = *fresh104;
            *fresh104 += 1;
            *SA
                .offset(
                    fresh105 as isize,
                ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        += 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(i as isize) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh106 = &mut (*induction_bucket
                .offset(
                    *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            let fresh107 = *fresh106;
            *fresh106 += 1;
            *SA
                .offset(
                    fresh107 as isize,
                ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_6k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let fresh108 = &mut (*buckets
        .offset(
            ((*T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 2 as std::ffi::c_int)
                + (*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    >= *T
                        .offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as fast_sint_t) as isize,
        ));
    let fresh109 = *fresh108;
    *fresh108 += 1;
    *SA
        .offset(
            fresh109 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | (-(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long);
    d += 1;
    *buckets
        .offset(
            (2 as std::ffi::c_int as std::ffi::c_long
                + ((*T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    << 2 as std::ffi::c_int)
                    + (*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                        >= *T
                            .offset(
                                (n - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            )) as std::ffi::c_int as fast_sint_t)) as isize,
        ) = d;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || left_suffixes_count < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        d = libsais16x64_partial_sorting_scan_left_to_right_32s_6k(
            T,
            SA,
            buckets,
            d,
            0 as std::ffi::c_int as fast_sint_t,
            left_suffixes_count,
        );
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_4k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset((2 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset((0 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let fresh110 = &mut (*induction_bucket
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh111 = *fresh110;
    *fresh110 += 1;
    *SA
        .offset(
            fresh111 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
            as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
        | (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
    d += 1;
    *distinct_names
        .offset(
            ((*T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as fast_sint_t) as isize,
        ) = d;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        d = libsais16x64_partial_sorting_scan_left_to_right_32s_4k(
            T,
            SA,
            k,
            buckets,
            d,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_left_to_right_32s_1k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let fresh112 = &mut (*buckets
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh113 = *fresh112;
    *fresh112 += 1;
    *SA
        .offset(
            fresh113 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
            as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_partial_sorting_scan_left_to_right_32s_1k(
            T,
            SA,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_shift_markers_16u_omp(
    mut SA: *mut sa_sint_t,
    mut _n: sa_sint_t,
    mut buckets: *const sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut temp_bucket: *const sa_sint_t = &*buckets
        .offset(
            (4 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *const sa_sint_t;
    let mut c: fast_sint_t = 0;
    c = (((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
        as fast_sint_t - 1 as std::ffi::c_int as std::ffi::c_long)
        << 1 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t;
    while c
        >= ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = 0;
        let mut s: sa_sint_t = -(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long;
        i = *temp_bucket.offset(c as isize) - 1 as std::ffi::c_int as std::ffi::c_long;
        j = *buckets
            .offset(
                (c
                    - (((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchw(
                &mut *SA.offset((i - prefetch_distance) as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            let mut p0: sa_sint_t = *SA
                .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q0: sa_sint_t = p0 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q0;
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p0 ^ q0;
            let mut p1: sa_sint_t = *SA
                .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q1: sa_sint_t = p1 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q1;
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p1 ^ q1;
            let mut p2: sa_sint_t = *SA
                .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q2: sa_sint_t = p2 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q2;
            *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p2 ^ q2;
            let mut p3: sa_sint_t = *SA
                .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q3: sa_sint_t = p3 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q3;
            *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p3 ^ q3;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            let mut p: sa_sint_t = *SA.offset(i as isize);
            let mut q: sa_sint_t = p & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q;
            *SA.offset(i as isize) = p ^ q;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c
            -= ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_shift_markers_32s_6k_omp(
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *const sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut temp_bucket: *const sa_sint_t = &*buckets
        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *const sa_sint_t;
    let mut c: fast_sint_t = 0;
    c = k - 1 as std::ffi::c_int as std::ffi::c_long;
    while c >= 1 as std::ffi::c_int as std::ffi::c_long {
        let mut i: fast_sint_t = 0;
        let mut j: fast_sint_t = 0;
        let mut s: sa_sint_t = -(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long;
        i = *buckets
            .offset(
                ((c << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t)
                    as isize,
            ) - 1 as std::ffi::c_int as std::ffi::c_long;
        j = *temp_bucket
            .offset(
                (((c - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t) as isize,
            ) + 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            libsais16x64_prefetchw(
                &mut *SA.offset((i - prefetch_distance) as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            let mut p0: sa_sint_t = *SA
                .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q0: sa_sint_t = p0 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q0;
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p0 ^ q0;
            let mut p1: sa_sint_t = *SA
                .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q1: sa_sint_t = p1 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q1;
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p1 ^ q1;
            let mut p2: sa_sint_t = *SA
                .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q2: sa_sint_t = p2 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q2;
            *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p2 ^ q2;
            let mut p3: sa_sint_t = *SA
                .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
            let mut q3: sa_sint_t = p3 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q3;
            *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = p3 ^ q3;
            i -= 4 as std::ffi::c_int as std::ffi::c_long;
        }
        j -= 3 as std::ffi::c_int as std::ffi::c_long;
        while i >= j {
            let mut p: sa_sint_t = *SA.offset(i as isize);
            let mut q: sa_sint_t = p & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) ^ s;
            s ^= q;
            *SA.offset(i as isize) = p ^ q;
            i -= 1 as std::ffi::c_int as std::ffi::c_long;
        }
        c -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_shift_markers_32s_4k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut s: sa_sint_t = (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
    i = n - 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchw(
            &mut *SA.offset((i - prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let mut q0: sa_sint_t = (p0
            & (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int)
            ^ s)
            & ((p0 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        s ^= q0;
        *SA.offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize) = p0 ^ q0;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        let mut q1: sa_sint_t = (p1
            & (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int)
            ^ s)
            & ((p1 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        s ^= q1;
        *SA.offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize) = p1 ^ q1;
        let mut p2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        let mut q2: sa_sint_t = (p2
            & (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int)
            ^ s)
            & ((p2 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        s ^= q2;
        *SA.offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize) = p2 ^ q2;
        let mut p3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        let mut q3: sa_sint_t = (p3
            & (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int)
            ^ s)
            & ((p3 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        s ^= q3;
        *SA.offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize) = p3 ^ q3;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        let mut q: sa_sint_t = (p
            & (1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int)
            ^ s)
            & ((p > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int);
        s ^= q;
        *SA.offset(i as isize) = p ^ q;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_shift_buckets_32s_6k(
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let mut temp_bucket: *mut sa_sint_t = &mut *buckets
        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    i = ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
        + 0 as std::ffi::c_int as fast_sint_t;
    while i
        <= ((k - 1 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t
    {
        *buckets
            .offset(
                (2 as std::ffi::c_int as std::ffi::c_long * i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = *temp_bucket
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        *buckets
            .offset(
                (2 as std::ffi::c_int as std::ffi::c_long * i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 2 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = *temp_bucket
            .offset(
                (i
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            );
        i
            += ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (0 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset(
            (2 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p0 &= 9223372036854775807 as std::ffi::c_long;
        let mut v0: sa_sint_t = ((*T
            .offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh114 = &mut (*induction_bucket.offset(v0 as isize));
        *fresh114 -= 1;
        *SA
            .offset(
                *fresh114 as isize,
            ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v0 as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v0 as isize) = d;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p1 &= 9223372036854775807 as std::ffi::c_long;
        let mut v1: sa_sint_t = ((*T
            .offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh115 = &mut (*induction_bucket.offset(v1 as isize));
        *fresh115 -= 1;
        *SA
            .offset(
                *fresh115 as isize,
            ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v1 as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v1 as isize) = d;
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        d
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p &= 9223372036854775807 as std::ffi::c_long;
        let mut v: sa_sint_t = ((*T
            .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        let fresh116 = &mut (*induction_bucket.offset(v as isize));
        *fresh116 -= 1;
        *SA
            .offset(
                *fresh116 as isize,
            ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v as isize) != d) as std::ffi::c_int
                as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_sint_t;
        *distinct_names.offset(v as isize) = d;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_gsa_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset(
            (0 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset(
            (2 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const uint16_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p0 &= 9223372036854775807 as std::ffi::c_long;
        let mut v0: sa_sint_t = ((*T
            .offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        if v0 != 1 as std::ffi::c_int as std::ffi::c_long {
            let fresh117 = &mut (*induction_bucket.offset(v0 as isize));
            *fresh117 -= 1;
            *SA
                .offset(
                    *fresh117 as isize,
                ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v0 as isize) != d) as std::ffi::c_int
                    as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                    as sa_sint_t;
            *distinct_names.offset(v0 as isize) = d;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p1 &= 9223372036854775807 as std::ffi::c_long;
        let mut v1: sa_sint_t = ((*T
            .offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        if v1 != 1 as std::ffi::c_int as std::ffi::c_long {
            let fresh118 = &mut (*induction_bucket.offset(v1 as isize));
            *fresh118 -= 1;
            *SA
                .offset(
                    *fresh118 as isize,
                ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v1 as isize) != d) as std::ffi::c_int
                    as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                    as sa_sint_t;
            *distinct_names.offset(v1 as isize) = d;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        d
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p &= 9223372036854775807 as std::ffi::c_long;
        let mut v: sa_sint_t = ((*T
            .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as fast_sint_t) << 1 as std::ffi::c_int)
            + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int
                > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as std::ffi::c_int) as std::ffi::c_int as fast_sint_t;
        if v != 1 as std::ffi::c_int as std::ffi::c_long {
            let fresh119 = &mut (*induction_bucket.offset(v as isize));
            *fresh119 -= 1;
            *SA
                .offset(
                    *fresh119 as isize,
                ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*distinct_names.offset(v as isize) != d) as std::ffi::c_int
                    as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int))
                    as sa_sint_t;
            *distinct_names.offset(v as isize) = d;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut scan_start: fast_sint_t = left_suffixes_count
        + 1 as std::ffi::c_int as std::ffi::c_long;
    let mut scan_end: fast_sint_t = n - first_lms_suffix;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || scan_end - scan_start < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_partial_sorting_scan_right_to_left_16u(
            T,
            SA,
            buckets,
            d,
            scan_start,
            scan_end - scan_start,
        );
    }
}
unsafe extern "C" fn libsais16x64_partial_gsa_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut scan_start: fast_sint_t = left_suffixes_count
        + 1 as std::ffi::c_int as std::ffi::c_long;
    let mut scan_end: fast_sint_t = n - first_lms_suffix;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || scan_end - scan_start < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_partial_gsa_scan_right_to_left_16u(
            T,
            SA,
            buckets,
            d,
            scan_start,
            scan_end - scan_start,
        );
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_6k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            (&*T
                .offset(
                    (*SA
                        .offset(
                            (i
                                - 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long) as isize,
                ) as *const sa_sint_t)
                .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) & 9223372036854775807 as std::ffi::c_long;
        let mut v0: sa_sint_t = (*T
            .offset(
                (p0
                    - (p0 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                        as std::ffi::c_long) as isize,
            ) << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t;
        libsais16x64_prefetchw(
            &mut *buckets.offset(v0 as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p1: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) & 9223372036854775807 as std::ffi::c_long;
        let mut v1: sa_sint_t = (*T
            .offset(
                (p1
                    - (p1 > 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                        as std::ffi::c_long) as isize,
            ) << 2 as std::ffi::c_int) + 0 as std::ffi::c_int as fast_sint_t;
        libsais16x64_prefetchw(
            &mut *buckets.offset(v1 as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p2: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p2 &= 9223372036854775807 as std::ffi::c_long;
        let mut v2: sa_sint_t = (*T
            .offset((p2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p2 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                > *T.offset((p2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh120 = &mut (*buckets.offset(v2 as isize));
        *fresh120 -= 1;
        *SA
            .offset(
                *fresh120 as isize,
            ) = (p2 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets
                .offset((2 as std::ffi::c_int as std::ffi::c_long + v2) as isize) != d)
                as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v2) as isize) = d;
        let mut p3: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        d
            += (p3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p3 &= 9223372036854775807 as std::ffi::c_long;
        let mut v3: sa_sint_t = (*T
            .offset((p3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p3 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                > *T.offset((p3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh121 = &mut (*buckets.offset(v3 as isize));
        *fresh121 -= 1;
        *SA
            .offset(
                *fresh121 as isize,
            ) = (p3 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets
                .offset((2 as std::ffi::c_int as std::ffi::c_long + v3) as isize) != d)
                as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v3) as isize) = d;
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        d
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p &= 9223372036854775807 as std::ffi::c_long;
        let mut v: sa_sint_t = (*T
            .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            << 2 as std::ffi::c_int)
            + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                as std::ffi::c_int as fast_sint_t;
        let fresh122 = &mut (*buckets.offset(v as isize));
        *fresh122 -= 1;
        *SA
            .offset(
                *fresh122 as isize,
            ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v) as isize)
                != d) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long + v) as isize) = d;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_4k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut induction_bucket: *mut sa_sint_t = &mut *buckets
        .offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut distinct_names: *mut sa_sint_t = &mut *buckets
        .offset((0 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int))
                } else {
                    2 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            Ts0.offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int))
                } else {
                    2 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            Ts1.offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            let Ts2: fast_sint_t = *T
                .offset(
                    ((s2
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int)))
                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            libsais16x64_prefetchw(
                &mut *induction_bucket.offset(Ts2 as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *distinct_names
                    .offset(
                        ((Ts2 << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            let Ts3: fast_sint_t = *T
                .offset(
                    ((s3
                        & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                                - 1 as std::ffi::c_int)))
                        - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            libsais16x64_prefetchw(
                &mut *induction_bucket.offset(Ts3 as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
            );
            libsais16x64_prefetchw(
                &mut *distinct_names
                    .offset(
                        ((Ts3 << 1 as std::ffi::c_int)
                            + 0 as std::ffi::c_int as fast_sint_t) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p0 >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p0
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v0: sa_sint_t = (*T
                .offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as fast_sint_t;
            let fresh123 = &mut (*induction_bucket
                .offset(
                    *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh123 -= 1;
            *SA
                .offset(
                    *fresh123 as isize,
                ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v0 as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v0 as isize) = d;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p1 >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p1
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v1: sa_sint_t = (*T
                .offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as fast_sint_t;
            let fresh124 = &mut (*induction_bucket
                .offset(
                    *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh124 -= 1;
            *SA
                .offset(
                    *fresh124 as isize,
                ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v1 as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v1 as isize) = d;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(i as isize) = 0 as std::ffi::c_int as sa_sint_t;
            d
                += p >> (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            p
                &= !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int));
            let mut v: sa_sint_t = (*T
                .offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                << 1 as std::ffi::c_int)
                + (*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as fast_sint_t;
            let fresh125 = &mut (*induction_bucket
                .offset(
                    *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh125 -= 1;
            *SA
                .offset(
                    *fresh125 as isize,
                ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t
                | ((*distinct_names.offset(v as isize) != d) as std::ffi::c_int
                    as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int);
            *distinct_names.offset(v as isize) = d;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_1k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s2 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s3 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh126 = &mut (*induction_bucket
                .offset(
                    *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh126 -= 1;
            *SA
                .offset(
                    *fresh126 as isize,
                ) = (p0 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p0 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh127 = &mut (*induction_bucket
                .offset(
                    *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh127 -= 1;
            *SA
                .offset(
                    *fresh127 as isize,
                ) = (p1 - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p1 - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T
                        .offset(
                            (p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        )) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(i as isize) = 0 as std::ffi::c_int as sa_sint_t;
            let fresh128 = &mut (*induction_bucket
                .offset(
                    *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                        as isize,
                ));
            *fresh128 -= 1;
            *SA
                .offset(
                    *fresh128 as isize,
                ) = (p - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((p - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    > *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_6k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut scan_start: fast_sint_t = left_suffixes_count
        + 1 as std::ffi::c_int as std::ffi::c_long;
    let mut scan_end: fast_sint_t = n - first_lms_suffix;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || scan_end - scan_start < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        d = libsais16x64_partial_sorting_scan_right_to_left_32s_6k(
            T,
            SA,
            buckets,
            d,
            scan_start,
            scan_end - scan_start,
        );
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_4k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut d: sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        d = libsais16x64_partial_sorting_scan_right_to_left_32s_4k(
            T,
            SA,
            k,
            buckets,
            d,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    d
}
unsafe extern "C" fn libsais16x64_partial_sorting_scan_right_to_left_32s_1k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_partial_sorting_scan_right_to_left_32s_1k(
            T,
            SA,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_gather_lms_suffixes_32s_4k(
    mut SA: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> fast_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut l: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - 3 as std::ffi::c_int as std::ffi::c_long;
    l = omp_block_start;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut s0: sa_uint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SA
            .offset(
                l as isize,
            ) = (s0
            .wrapping_sub(
                ((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int)) as sa_uint_t,
            )
            & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_uint_t) as sa_sint_t;
        l
            += ((s0 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut s1: sa_uint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SA
            .offset(
                l as isize,
            ) = (s1
            .wrapping_sub(
                ((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int)) as sa_uint_t,
            )
            & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_uint_t) as sa_sint_t;
        l
            += ((s1 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut s2: sa_uint_t = *SA
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SA
            .offset(
                l as isize,
            ) = (s2
            .wrapping_sub(
                ((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int)) as sa_uint_t,
            )
            & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_uint_t) as sa_sint_t;
        l
            += ((s2 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut s3: sa_uint_t = *SA
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SA
            .offset(
                l as isize,
            ) = (s3
            .wrapping_sub(
                ((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int)) as sa_uint_t,
            )
            & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_uint_t) as sa_sint_t;
        l
            += ((s3 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut s: sa_uint_t = *SA.offset(i as isize) as sa_uint_t;
        *SA
            .offset(
                l as isize,
            ) = (s
            .wrapping_sub(
                ((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int
                        - 1 as std::ffi::c_int)) as sa_uint_t,
            )
            & !((1 as std::ffi::c_int as sa_sint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int - 1 as std::ffi::c_int))
                as sa_uint_t) as sa_sint_t;
        l
            += ((s as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    l
}
unsafe extern "C" fn libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k(
    mut SA: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> fast_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut l: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - 3 as std::ffi::c_int as std::ffi::c_long;
    l = omp_block_start;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s0 & 9223372036854775807 as std::ffi::c_long;
        l
            += (s0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s1 & 9223372036854775807 as std::ffi::c_long;
        l
            += (s1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s2: sa_sint_t = *SA
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s2 & 9223372036854775807 as std::ffi::c_long;
        l
            += (s2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s3: sa_sint_t = *SA
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s3 & 9223372036854775807 as std::ffi::c_long;
        l
            += (s3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut s: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(l as isize) = s & 9223372036854775807 as std::ffi::c_long;
        l
            += (s < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    l
}
unsafe extern "C" fn libsais16x64_partial_sorting_gather_lms_suffixes_32s_4k_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_partial_sorting_gather_lms_suffixes_32s_4k(
            SA,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k(
            SA,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_induce_partial_order_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut flags: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    memset(
        &mut *buckets
            .offset(
                (2 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as isize,
            ) as *mut sa_sint_t as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (2 as std::ffi::c_int as size_t)
            .wrapping_mul(
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong,
            )
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
        *buckets
            .offset(
                ((4 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as std::ffi::c_long
                    + (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) = *buckets
            .offset(
                ((4 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as std::ffi::c_long
                    + (((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 1 as std::ffi::c_int as fast_sint_t)) as isize,
            ) - 1 as std::ffi::c_int as std::ffi::c_long;
        libsais16x64_flip_suffix_markers_omp(
            SA,
            *buckets
                .offset(
                    ((4 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as std::ffi::c_long
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                ),
            threads,
        );
    }
    let mut d: sa_sint_t = libsais16x64_partial_sorting_scan_left_to_right_16u_omp(
        T,
        SA,
        n,
        k,
        buckets,
        left_suffixes_count,
        0 as std::ffi::c_int as sa_sint_t,
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_shift_markers_16u_omp(SA, n, buckets, threads);
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
        libsais16x64_partial_gsa_scan_right_to_left_16u_omp(
            T,
            SA,
            n,
            k,
            buckets,
            first_lms_suffix,
            left_suffixes_count,
            d,
            threads,
            thread_state,
        );
        if *T.offset(first_lms_suffix as isize) as std::ffi::c_int
            == 0 as std::ffi::c_int
        {
            memmove(
                &mut *SA.offset(1 as std::ffi::c_int as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
                ((*buckets
                    .offset(
                        (((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t) as isize,
                    ) - 1 as std::ffi::c_int as std::ffi::c_long) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            *SA
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = first_lms_suffix | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        }
        *buckets
            .offset(
                (((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                    + 1 as std::ffi::c_int as fast_sint_t) as isize,
            ) = 0 as std::ffi::c_int as sa_sint_t;
    } else {
        libsais16x64_partial_sorting_scan_right_to_left_16u_omp(
            T,
            SA,
            n,
            k,
            buckets,
            first_lms_suffix,
            left_suffixes_count,
            d,
            threads,
            thread_state,
        );
    };
}
unsafe extern "C" fn libsais16x64_induce_partial_order_32s_6k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut first_lms_suffix: sa_sint_t,
    mut left_suffixes_count: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut d: sa_sint_t = libsais16x64_partial_sorting_scan_left_to_right_32s_6k_omp(
        T,
        SA,
        n,
        buckets,
        left_suffixes_count,
        0 as std::ffi::c_int as sa_sint_t,
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_shift_markers_32s_6k_omp(SA, k, buckets, threads);
    libsais16x64_partial_sorting_shift_buckets_32s_6k(k, buckets);
    libsais16x64_partial_sorting_scan_right_to_left_32s_6k_omp(
        T,
        SA,
        n,
        buckets,
        first_lms_suffix,
        left_suffixes_count,
        d,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_partial_order_32s_4k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    memset(
        buckets as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (2 as std::ffi::c_int as std::ffi::c_ulong)
            .wrapping_mul(k as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut d: sa_sint_t = libsais16x64_partial_sorting_scan_left_to_right_32s_4k_omp(
        T,
        SA,
        n,
        k,
        buckets,
        0 as std::ffi::c_int as sa_sint_t,
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_shift_markers_32s_4k(SA, n);
    libsais16x64_partial_sorting_scan_right_to_left_32s_4k_omp(
        T,
        SA,
        n,
        k,
        buckets,
        d,
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_gather_lms_suffixes_32s_4k_omp(
        SA,
        n,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_partial_order_32s_2k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_partial_sorting_scan_left_to_right_32s_1k_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((1 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_scan_right_to_left_32s_1k_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((0 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k_omp(
        SA,
        n,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_partial_order_32s_1k_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_start_32s_1k(k, buckets);
    libsais16x64_partial_sorting_scan_left_to_right_32s_1k_omp(
        T,
        SA,
        n,
        buckets,
        threads,
        thread_state,
    );
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_end_32s_1k(k, buckets);
    libsais16x64_partial_sorting_scan_right_to_left_32s_1k_omp(
        T,
        SA,
        n,
        buckets,
        threads,
        thread_state,
    );
    libsais16x64_partial_sorting_gather_lms_suffixes_32s_1k_omp(
        SA,
        n,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_renumber_lms_suffixes_16u(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut name: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAm: *mut sa_sint_t = &mut *SA.offset(m as isize) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SAm
            .offset(
                ((p0 & 9223372036854775807 as std::ffi::c_long) >> 1 as std::ffi::c_int)
                    as isize,
            ) = name | (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SAm
            .offset(
                ((p1 & 9223372036854775807 as std::ffi::c_long) >> 1 as std::ffi::c_int)
                    as isize,
            ) = name | (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut p2: sa_sint_t = *SA
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SAm
            .offset(
                ((p2 & 9223372036854775807 as std::ffi::c_long) >> 1 as std::ffi::c_int)
                    as isize,
            ) = name | (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut p3: sa_sint_t = *SA
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SAm
            .offset(
                ((p3 & 9223372036854775807 as std::ffi::c_long) >> 1 as std::ffi::c_int)
                    as isize,
            ) = name | (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SAm
            .offset(
                ((p & 9223372036854775807 as std::ffi::c_long) >> 1 as std::ffi::c_int)
                    as isize,
            ) = name | (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    name
}
unsafe extern "C" fn libsais16x64_gather_marked_lms_suffixes(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut l: fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> fast_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    l -= 1 as std::ffi::c_int as std::ffi::c_long;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = m + omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = m + omp_block_start + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA.offset((i - prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s0 & 9223372036854775807 as std::ffi::c_long;
        l
            -= (s0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s1 & 9223372036854775807 as std::ffi::c_long;
        l
            -= (s1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s2 & 9223372036854775807 as std::ffi::c_long;
        l
            -= (s2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        let mut s3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA.offset(l as isize) = s3 & 9223372036854775807 as std::ffi::c_long;
        l
            -= (s3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut s: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(l as isize) = s & 9223372036854775807 as std::ffi::c_long;
        l
            -= (s < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    l += 1 as std::ffi::c_int as std::ffi::c_long;
    l
}
unsafe extern "C" fn libsais16x64_renumber_lms_suffixes_16u_omp(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut name: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (m / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        m - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        name = libsais16x64_renumber_lms_suffixes_16u(
            SA,
            m,
            0 as std::ffi::c_int as sa_sint_t,
            omp_block_start,
            omp_block_size,
        );
    }
    name
}
unsafe extern "C" fn libsais16x64_gather_marked_lms_suffixes_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = ((n >> 1 as std::ffi::c_int) / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        (n >> 1 as std::ffi::c_int) - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_gather_marked_lms_suffixes(
            SA,
            m,
            n + fs,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_renumber_and_gather_lms_suffixes_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    memset(
        &mut *SA.offset(m as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (n as size_t >> 1 as std::ffi::c_int)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut name: sa_sint_t = libsais16x64_renumber_lms_suffixes_16u_omp(
        SA,
        m,
        threads,
        thread_state,
    );
    if name < m {
        libsais16x64_gather_marked_lms_suffixes_omp(SA, n, m, fs, threads, thread_state);
    } else {
        let mut i: fast_sint_t = 0;
        i = 0 as std::ffi::c_int as fast_sint_t;
        while i < m {
            let fresh129 = &mut (*SA.offset(i as isize));
            *fresh129 &= 9223372036854775807 as std::ffi::c_long;
            i += 1 as std::ffi::c_int as std::ffi::c_long;
        }
    }
    name
}
unsafe extern "C" fn libsais16x64_renumber_distinct_lms_suffixes_32s_4k(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut name: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAm: *mut sa_sint_t = &mut *SA.offset(m as isize) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut p0: sa_sint_t = 0;
    let mut p1: sa_sint_t = 0;
    let mut p2: sa_sint_t = 0;
    let mut p3: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    ((*SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) & 9223372036854775807 as std::ffi::c_long)
                        >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        p0 = *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh130 = &mut (*SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize));
        *fresh130 = p0 & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (*fresh130 >> 1 as std::ffi::c_int) as isize,
            ) = name
            | p0 & p3 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p1 = *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh131 = &mut (*SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize));
        *fresh131 = p1 & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (*fresh131 >> 1 as std::ffi::c_int) as isize,
            ) = name
            | p1 & p0 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p1 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p2 = *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh132 = &mut (*SA
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize));
        *fresh132 = p2 & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (*fresh132 >> 1 as std::ffi::c_int) as isize,
            ) = name
            | p2 & p1 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p2 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p3 = *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        let fresh133 = &mut (*SA
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize));
        *fresh133 = p3 & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (*fresh133 >> 1 as std::ffi::c_int) as isize,
            ) = name
            | p3 & p2 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        p2 = p3;
        p3 = *SA.offset(i as isize);
        let fresh134 = &mut (*SA.offset(i as isize));
        *fresh134 = p3 & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (*fresh134 >> 1 as std::ffi::c_int) as isize,
            ) = name
            | p3 & p2 & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        name
            += (p3 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    name
}
unsafe extern "C" fn libsais16x64_mark_distinct_lms_suffixes_32s(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut p0: sa_sint_t = 0;
    let mut p1: sa_sint_t = 0;
    let mut p2: sa_sint_t = 0;
    let mut p3: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    i = m + omp_block_start;
    j = m + omp_block_start + omp_block_size - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        p0 = *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & (p3 | 9223372036854775807 as std::ffi::c_long);
        p0 = if p0 == 0 as std::ffi::c_int as std::ffi::c_long { p3 } else { p0 };
        p1 = *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & (p0 | 9223372036854775807 as std::ffi::c_long);
        p1 = if p1 == 0 as std::ffi::c_int as std::ffi::c_long { p0 } else { p1 };
        p2 = *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p2 & (p1 | 9223372036854775807 as std::ffi::c_long);
        p2 = if p2 == 0 as std::ffi::c_int as std::ffi::c_long { p1 } else { p2 };
        p3 = *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p3 & (p2 | 9223372036854775807 as std::ffi::c_long);
        p3 = if p3 == 0 as std::ffi::c_int as std::ffi::c_long { p2 } else { p3 };
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        p2 = p3;
        p3 = *SA.offset(i as isize);
        *SA.offset(i as isize) = p3 & (p2 | 9223372036854775807 as std::ffi::c_long);
        p3 = if p3 == 0 as std::ffi::c_int as std::ffi::c_long { p2 } else { p3 };
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_clamp_lms_suffixes_length_32s(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAm: *mut sa_sint_t = &mut *SA.offset(m as isize) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SAm.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        *SAm
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = (if *SAm.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            < 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SAm.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
        } else {
            0 as std::ffi::c_int as std::ffi::c_long
        }) & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = (if *SAm.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            < 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SAm.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
        } else {
            0 as std::ffi::c_int as std::ffi::c_long
        }) & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = (if *SAm.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            < 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SAm.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
        } else {
            0 as std::ffi::c_int as std::ffi::c_long
        }) & 9223372036854775807 as std::ffi::c_long;
        *SAm
            .offset(
                (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = (if *SAm.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            < 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SAm.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
        } else {
            0 as std::ffi::c_int as std::ffi::c_long
        }) & 9223372036854775807 as std::ffi::c_long;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *SAm
            .offset(
                i as isize,
            ) = (if *SAm.offset(i as isize) < 0 as std::ffi::c_int as std::ffi::c_long {
            *SAm.offset(i as isize)
        } else {
            0 as std::ffi::c_int as std::ffi::c_long
        }) & 9223372036854775807 as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_renumber_distinct_lms_suffixes_32s_4k_omp(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut name: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (m / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        m - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        name = libsais16x64_renumber_distinct_lms_suffixes_32s_4k(
            SA,
            m,
            1 as std::ffi::c_int as sa_sint_t,
            omp_block_start,
            omp_block_size,
        );
    }
    name - 1 as std::ffi::c_int as std::ffi::c_long
}
unsafe extern "C" fn libsais16x64_mark_distinct_lms_suffixes_32s_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = n >> 1 as std::ffi::c_int;
    libsais16x64_mark_distinct_lms_suffixes_32s(SA, m, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_clamp_lms_suffixes_length_32s_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = n >> 1 as std::ffi::c_int;
    libsais16x64_clamp_lms_suffixes_length_32s(SA, m, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    memset(
        &mut *SA.offset(m as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (n as size_t >> 1 as std::ffi::c_int)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut name: sa_sint_t = libsais16x64_renumber_distinct_lms_suffixes_32s_4k_omp(
        SA,
        m,
        threads,
        thread_state,
    );
    if name < m {
        libsais16x64_mark_distinct_lms_suffixes_32s_omp(SA, n, m, threads);
    }
    name
}
unsafe extern "C" fn libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut threads: sa_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAm: *mut sa_sint_t = &mut *SA.offset(m as isize) as *mut sa_sint_t;
    libsais16x64_gather_lms_suffixes_32s(T, SA, n);
    memset(
        &mut *SA.offset(m as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (n as size_t)
            .wrapping_sub(m as size_t)
            .wrapping_sub(m as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = n - m;
    j = n - 1 as std::ffi::c_int as std::ffi::c_long - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        *SAm
            .offset(
                (*SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as sa_uint_t >> 1 as std::ffi::c_int) as isize,
            ) = *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            - *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            + 1 as std::ffi::c_int as std::ffi::c_long
            + (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        *SAm
            .offset(
                (*SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as sa_uint_t >> 1 as std::ffi::c_int) as isize,
            ) = *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            - *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            + 1 as std::ffi::c_int as std::ffi::c_long
            + (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        *SAm
            .offset(
                (*SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as sa_uint_t >> 1 as std::ffi::c_int) as isize,
            ) = *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            - *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            + 1 as std::ffi::c_int as std::ffi::c_long
            + (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        *SAm
            .offset(
                (*SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as sa_uint_t >> 1 as std::ffi::c_int) as isize,
            ) = *SA.offset((i + 4 as std::ffi::c_int as std::ffi::c_long) as isize)
            - *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            + 1 as std::ffi::c_int as std::ffi::c_long
            + (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *SAm
            .offset(
                (*SA.offset(i as isize) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
            ) = *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            - *SA.offset(i as isize) + 1 as std::ffi::c_int as std::ffi::c_long
            + (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    *SAm
        .offset(
            (*SA.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as sa_uint_t >> 1 as std::ffi::c_int) as isize,
        ) = 1 as std::ffi::c_int as std::ffi::c_long
        + (-(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long);
    libsais16x64_clamp_lms_suffixes_length_32s_omp(SA, n, m, threads);
    let mut name: sa_sint_t = 1 as std::ffi::c_int as sa_sint_t;
    let mut i_0: fast_sint_t = 0;
    let mut j_0: fast_sint_t = 0;
    let mut p: fast_sint_t = *SA.offset(0 as std::ffi::c_int as isize);
    let mut plen: fast_sint_t = *SAm.offset((p >> 1 as std::ffi::c_int) as isize);
    let mut pdiff: sa_sint_t = -(9223372036854775807 as std::ffi::c_long)
        - 1 as std::ffi::c_int as std::ffi::c_long;
    i_0 = 1 as std::ffi::c_int as fast_sint_t;
    j_0 = m - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long;
    while i_0 < j_0 {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i_0 + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i_0 + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &mut *T
                .offset(
                    *SA
                        .offset(
                            (i_0 + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i_0 + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &mut *T
                .offset(
                    *SA
                        .offset(
                            (i_0 + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut q: fast_sint_t = *SA
            .offset((i_0 + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        let mut qlen: fast_sint_t = *SAm.offset((q >> 1 as std::ffi::c_int) as isize);
        let mut qdiff: sa_sint_t = -(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long;
        if plen == qlen {
            let mut l: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
            while *T.offset((p + l) as isize) == *T.offset((q + l) as isize) {
                l += 1;
                if l >= qlen {
                    break;
                }
            }
            qdiff = (l - qlen) & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        }
        *SAm.offset((p >> 1 as std::ffi::c_int) as isize) = name | pdiff & qdiff;
        name
            += (qdiff < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p = *SA.offset((i_0 + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        plen = *SAm.offset((p >> 1 as std::ffi::c_int) as isize);
        pdiff = -(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long;
        if qlen == plen {
            let mut l_0: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
            while *T.offset((q + l_0) as isize) == *T.offset((p + l_0) as isize) {
                l_0 += 1;
                if l_0 >= plen {
                    break;
                }
            }
            pdiff = (l_0 - plen) & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        }
        *SAm.offset((q >> 1 as std::ffi::c_int) as isize) = name | qdiff & pdiff;
        name
            += (pdiff < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i_0 += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j_0 += prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i_0 < j_0 {
        let mut q_0: fast_sint_t = *SA.offset(i_0 as isize);
        let mut qlen_0: fast_sint_t = *SAm
            .offset((q_0 >> 1 as std::ffi::c_int) as isize);
        let mut qdiff_0: sa_sint_t = -(9223372036854775807 as std::ffi::c_long)
            - 1 as std::ffi::c_int as std::ffi::c_long;
        if plen == qlen_0 {
            let mut l_1: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
            while *T.offset((p + l_1) as isize) == *T.offset((q_0 + l_1) as isize) {
                l_1 += 1;
                if l_1 >= plen {
                    break;
                }
            }
            qdiff_0 = (l_1 - plen) & (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
        }
        *SAm.offset((p >> 1 as std::ffi::c_int) as isize) = name | pdiff & qdiff_0;
        name
            += (qdiff_0 < 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        p = q_0;
        plen = qlen_0;
        pdiff = qdiff_0;
        i_0 += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    *SAm.offset((p >> 1 as std::ffi::c_int) as isize) = name | pdiff;
    name += 1;
    if name <= m {
        libsais16x64_mark_distinct_lms_suffixes_32s_omp(SA, n, m, threads);
    }
    name - 1 as std::ffi::c_int as std::ffi::c_long
}
unsafe extern "C" fn libsais16x64_reconstruct_lms_suffixes(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAnm: *const sa_sint_t = &mut *SA.offset((n - m) as isize) as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*SAnm
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*SAnm
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*SAnm
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*SAnm
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *SAnm
            .offset(
                *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *SAnm
            .offset(
                *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *SAnm
            .offset(
                *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        *SA
            .offset(
                (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *SAnm
            .offset(
                *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *SA.offset(i as isize) = *SAnm.offset(*SA.offset(i as isize) as isize);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_reconstruct_lms_suffixes_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = m;
    libsais16x64_reconstruct_lms_suffixes(SA, n, m, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_interval_16u(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut flags: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
        let fresh135 = &mut (*buckets
            .offset(
                (7 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as isize,
            ));
        *fresh135 -= 1;
    }
    let mut bucket_end: *const sa_sint_t = &mut *buckets
        .offset(
            (7 as std::ffi::c_int
                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int)) as isize,
        ) as *mut sa_sint_t;
    let mut c: fast_sint_t = 0;
    let mut j: fast_sint_t = n;
    c = ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
        - 2 as std::ffi::c_int) as fast_sint_t;
    while c >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut l: fast_sint_t = *buckets
            .offset(
                ((c << 1 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t
                    + (((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            )
            - *buckets
                .offset(
                    ((c << 1 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t)
                        as isize,
                );
        if l > 0 as std::ffi::c_int as std::ffi::c_long {
            let mut i: fast_sint_t = *bucket_end.offset(c as isize);
            if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                memset(
                    &mut *SA.offset(i as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    0 as std::ffi::c_int,
                    ((j - i) as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            j = i - l;
            m -= l;
            memmove(
                &mut *SA.offset(j as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
                &mut *SA.offset(m as isize) as *mut sa_sint_t as *const std::ffi::c_void,
                (l as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        c -= 1;
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
        let fresh136 = &mut (*buckets
            .offset(
                (7 as std::ffi::c_int
                    * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int)) as isize,
            ));
        *fresh136 += 1;
    }
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_interval_32s_4k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *const sa_sint_t,
) {
    let mut bucket_end: *const sa_sint_t = &*buckets
        .offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *const sa_sint_t;
    let mut c: fast_sint_t = 0;
    let mut j: fast_sint_t = n;
    c = k - 2 as std::ffi::c_int as std::ffi::c_long;
    while c >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut l: fast_sint_t = *buckets
            .offset(
                ((c << 1 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t
                    + (((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                        + 0 as std::ffi::c_int as fast_sint_t)) as isize,
            )
            - *buckets
                .offset(
                    ((c << 1 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t)
                        as isize,
                );
        if l > 0 as std::ffi::c_int as std::ffi::c_long {
            let mut i: fast_sint_t = *bucket_end.offset(c as isize);
            if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                memset(
                    &mut *SA.offset(i as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    0 as std::ffi::c_int,
                    ((j - i) as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            j = i - l;
            m -= l;
            memmove(
                &mut *SA.offset(j as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
                &mut *SA.offset(m as isize) as *mut sa_sint_t as *const std::ffi::c_void,
                (l as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        c -= 1;
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_interval_32s_2k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *const sa_sint_t,
) {
    let mut j: fast_sint_t = n;
    if k > 1 as std::ffi::c_int as std::ffi::c_long {
        let mut c: fast_sint_t = 0;
        c = ((k - 2 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t;
        while c
            >= ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t
        {
            let mut l: fast_sint_t = *buckets
                .offset(
                    (c
                        + (((1 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                )
                - *buckets
                    .offset(
                        (c
                            + (((0 as std::ffi::c_int as fast_sint_t)
                                << 1 as std::ffi::c_int)
                                + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                    );
            if l > 0 as std::ffi::c_int as std::ffi::c_long {
                let mut i: fast_sint_t = *buckets.offset(c as isize);
                if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                    memset(
                        &mut *SA.offset(i as isize) as *mut sa_sint_t
                            as *mut std::ffi::c_void,
                        0 as std::ffi::c_int,
                        ((j - i) as size_t)
                            .wrapping_mul(
                                ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                            ),
                    );
                }
                j = i - l;
                m -= l;
                memmove(
                    &mut *SA.offset(j as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    &mut *SA.offset(m as isize) as *mut sa_sint_t
                        as *const std::ffi::c_void,
                    (l as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            c
                -= ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t;
        }
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_interval_32s_1k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *mut sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut c: sa_sint_t = k - 1 as std::ffi::c_int as std::ffi::c_long;
    let mut i: fast_sint_t = 0;
    let mut l: fast_sint_t = *buckets.offset(c as isize);
    i = m - 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T
                .offset(
                    *SA
                        .offset(
                            (i - prefetch_distance
                                - 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        if *T.offset(p0 as isize) != c {
            c = *T.offset(p0 as isize);
            memset(
                &mut *SA.offset(*buckets.offset(c as isize) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((l - *buckets.offset(c as isize)) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            l = *buckets.offset(c as isize);
        }
        l -= 1;
        *SA.offset(l as isize) = p0;
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        if *T.offset(p1 as isize) != c {
            c = *T.offset(p1 as isize);
            memset(
                &mut *SA.offset(*buckets.offset(c as isize) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((l - *buckets.offset(c as isize)) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            l = *buckets.offset(c as isize);
        }
        l -= 1;
        *SA.offset(l as isize) = p1;
        let mut p2: sa_sint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        if *T.offset(p2 as isize) != c {
            c = *T.offset(p2 as isize);
            memset(
                &mut *SA.offset(*buckets.offset(c as isize) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((l - *buckets.offset(c as isize)) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            l = *buckets.offset(c as isize);
        }
        l -= 1;
        *SA.offset(l as isize) = p2;
        let mut p3: sa_sint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        if *T.offset(p3 as isize) != c {
            c = *T.offset(p3 as isize);
            memset(
                &mut *SA.offset(*buckets.offset(c as isize) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((l - *buckets.offset(c as isize)) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            l = *buckets.offset(c as isize);
        }
        l -= 1;
        *SA.offset(l as isize) = p3;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    while i >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        if *T.offset(p as isize) != c {
            c = *T.offset(p as isize);
            memset(
                &mut *SA.offset(*buckets.offset(c as isize) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((l - *buckets.offset(c as isize)) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            l = *buckets.offset(c as isize);
        }
        l -= 1;
        *SA.offset(l as isize) = p;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (l as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_histogram_32s_6k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *const sa_sint_t,
) {
    let mut bucket_end: *const sa_sint_t = &*buckets
        .offset((5 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *const sa_sint_t;
    let mut c: fast_sint_t = 0;
    let mut j: fast_sint_t = n;
    c = k - 2 as std::ffi::c_int as std::ffi::c_long;
    while c >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut l: fast_sint_t = *buckets
            .offset(
                ((c << 2 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t)
                    as isize,
            );
        if l > 0 as std::ffi::c_int as std::ffi::c_long {
            let mut i: fast_sint_t = *bucket_end.offset(c as isize);
            if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                memset(
                    &mut *SA.offset(i as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    0 as std::ffi::c_int,
                    ((j - i) as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            j = i - l;
            m -= l;
            memmove(
                &mut *SA.offset(j as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
                &mut *SA.offset(m as isize) as *mut sa_sint_t as *const std::ffi::c_void,
                (l as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        c -= 1;
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_histogram_32s_4k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *const sa_sint_t,
) {
    let mut bucket_end: *const sa_sint_t = &*buckets
        .offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize)
        as *const sa_sint_t;
    let mut c: fast_sint_t = 0;
    let mut j: fast_sint_t = n;
    c = k - 2 as std::ffi::c_int as std::ffi::c_long;
    while c >= 0 as std::ffi::c_int as std::ffi::c_long {
        let mut l: fast_sint_t = *buckets
            .offset(
                ((c << 1 as std::ffi::c_int) + 1 as std::ffi::c_int as fast_sint_t)
                    as isize,
            );
        if l > 0 as std::ffi::c_int as std::ffi::c_long {
            let mut i: fast_sint_t = *bucket_end.offset(c as isize);
            if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                memset(
                    &mut *SA.offset(i as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    0 as std::ffi::c_int,
                    ((j - i) as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            j = i - l;
            m -= l;
            memmove(
                &mut *SA.offset(j as isize) as *mut sa_sint_t as *mut std::ffi::c_void,
                &mut *SA.offset(m as isize) as *mut sa_sint_t as *const std::ffi::c_void,
                (l as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        c -= 1;
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_place_lms_suffixes_histogram_32s_2k(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut buckets: *const sa_sint_t,
) {
    let mut j: fast_sint_t = n;
    if k > 1 as std::ffi::c_int as std::ffi::c_long {
        let mut c: fast_sint_t = 0;
        c = ((k - 2 as std::ffi::c_int as std::ffi::c_long) << 1 as std::ffi::c_int)
            + 0 as std::ffi::c_int as fast_sint_t;
        while c
            >= ((0 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                + 0 as std::ffi::c_int as fast_sint_t
        {
            let mut l: fast_sint_t = *buckets
                .offset(
                    (c
                        + (((0 as std::ffi::c_int as fast_sint_t)
                            << 1 as std::ffi::c_int)
                            + 1 as std::ffi::c_int as fast_sint_t)) as isize,
                );
            if l > 0 as std::ffi::c_int as std::ffi::c_long {
                let mut i: fast_sint_t = *buckets.offset(c as isize);
                if j - i > 0 as std::ffi::c_int as std::ffi::c_long {
                    memset(
                        &mut *SA.offset(i as isize) as *mut sa_sint_t
                            as *mut std::ffi::c_void,
                        0 as std::ffi::c_int,
                        ((j - i) as size_t)
                            .wrapping_mul(
                                ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                            ),
                    );
                }
                j = i - l;
                m -= l;
                memmove(
                    &mut *SA.offset(j as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    &mut *SA.offset(m as isize) as *mut sa_sint_t
                        as *const std::ffi::c_void,
                    (l as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            c
                -= ((1 as std::ffi::c_int as fast_sint_t) << 1 as std::ffi::c_int)
                    + 0 as std::ffi::c_int as fast_sint_t;
        }
    }
    memset(
        &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        0 as std::ffi::c_int,
        (j as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_final_bwt_scan_left_to_right_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            *SA
                .offset(
                    (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = *T.offset(p0 as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh137 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            let fresh138 = *fresh137;
            *fresh137 += 1;
            *SA
                .offset(
                    fresh138 as isize,
                ) = p0
                | ((((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p0 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            *SA
                .offset(
                    (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = *T.offset(p1 as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh139 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            let fresh140 = *fresh139;
            *fresh139 += 1;
            *SA
                .offset(
                    fresh140 as isize,
                ) = p1
                | ((((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p1 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            *SA
                .offset(
                    i as isize,
                ) = *T.offset(p as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh141 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            let fresh142 = *fresh141;
            *fresh141 += 1;
            *SA
                .offset(
                    fresh142 as isize,
                ) = p
                | ((((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_bwt_aux_scan_left_to_right_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut rm: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            *SA
                .offset(
                    (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = *T.offset(p0 as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh143 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            let fresh144 = *fresh143;
            *fresh143 += 1;
            *SA
                .offset(
                    fresh144 as isize,
                ) = p0
                | ((((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p0 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
            if p0 & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p0 / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p0 as isize) as isize);
            }
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            *SA
                .offset(
                    (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = *T.offset(p1 as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh145 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            let fresh146 = *fresh145;
            *fresh145 += 1;
            *SA
                .offset(
                    fresh146 as isize,
                ) = p1
                | ((((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p1 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
            if p1 & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p1 / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p1 as isize) as isize);
            }
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            *SA
                .offset(
                    i as isize,
                ) = *T.offset(p as isize) as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh147 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            let fresh148 = *fresh147;
            *fresh147 += 1;
            *SA
                .offset(
                    fresh148 as isize,
                ) = p
                | ((((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
            if p & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p as isize) as isize);
            }
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_left_to_right_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let fresh149 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            let fresh150 = *fresh149;
            *fresh149 += 1;
            *SA
                .offset(
                    fresh150 as isize,
                ) = p0
                | ((((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p0 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let fresh151 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            let fresh152 = *fresh151;
            *fresh151 += 1;
            *SA
                .offset(
                    fresh152 as isize,
                ) = p1
                | ((((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p1 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA
            .offset(
                i as isize,
            ) = p ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let fresh153 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            let fresh154 = *fresh153;
            *fresh153 += 1;
            *SA
                .offset(
                    fresh154 as isize,
                ) = p
                | ((((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int) < *T.offset(p as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_left_to_right_32s(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size
        - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        - 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s2 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s3 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let fresh155 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            let fresh156 = *fresh155;
            *fresh155 += 1;
            *SA
                .offset(
                    fresh156 as isize,
                ) = p0
                | (((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) < *T.offset(p0 as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let fresh157 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            let fresh158 = *fresh157;
            *fresh157 += 1;
            *SA
                .offset(
                    fresh158 as isize,
                ) = p1
                | (((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) < *T.offset(p1 as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        += 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA
            .offset(
                i as isize,
            ) = p ^ (-(9223372036854775807 as std::ffi::c_long)
                - 1 as std::ffi::c_int as std::ffi::c_long);
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let fresh159 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            let fresh160 = *fresh159;
            *fresh159 += 1;
            *SA
                .offset(
                    fresh160 as isize,
                ) = p
                | (((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) < *T.offset(p as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_bwt_scan_left_to_right_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: fast_sint_t,
    mut _k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let fresh161 = &mut (*induction_bucket
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh162 = *fresh161;
    *fresh161 += 1;
    *SA
        .offset(
            fresh162 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | ((((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as std::ffi::c_int)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_bwt_scan_left_to_right_16u(
            T,
            SA,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_bwt_aux_scan_left_to_right_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: fast_sint_t,
    mut _k: sa_sint_t,
    mut rm: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let fresh163 = &mut (*induction_bucket
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh164 = *fresh163;
    *fresh163 += 1;
    *SA
        .offset(
            fresh164 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | ((((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as std::ffi::c_int)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
    if (n - 1 as std::ffi::c_int as std::ffi::c_long) & rm
        == 0 as std::ffi::c_int as std::ffi::c_long
    {
        *I
            .offset(
                ((n - 1 as std::ffi::c_int as std::ffi::c_long)
                    / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
            ) = *induction_bucket
            .offset(
                *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
    }
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_bwt_aux_scan_left_to_right_16u(
            T,
            SA,
            rm,
            I,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_left_to_right_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: fast_sint_t,
    mut _k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let fresh165 = &mut (*induction_bucket
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh166 = *fresh165;
    *fresh165 += 1;
    *SA
        .offset(
            fresh166 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | ((((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as std::ffi::c_int)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_sorting_scan_left_to_right_16u(
            T,
            SA,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_left_to_right_32s_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let fresh167 = &mut (*induction_bucket
        .offset(
            *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize) as isize,
        ));
    let fresh168 = *fresh167;
    *fresh167 += 1;
    *SA
        .offset(
            fresh168 as isize,
        ) = (n - 1 as std::ffi::c_int as std::ffi::c_long) | (((*T.offset((n - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            < *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize))
            as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_sorting_scan_left_to_right_32s(
            T,
            SA,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_bwt_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut index: sa_sint_t = -(1 as std::ffi::c_int) as sa_sint_t;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        index = if p0 == 0 as std::ffi::c_int as std::ffi::c_long {
            i - 0 as std::ffi::c_int as std::ffi::c_long
        } else {
            index
        };
        *SA
            .offset(
                (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let mut c0: uint16_t = *T
                .offset(
                    (p0
                        - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1: uint16_t = *T.offset(p0 as isize);
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c1 as sa_sint_t;
            let mut t: sa_sint_t = c0 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh169 = &mut (*induction_bucket.offset(c1 as isize));
            *fresh169 -= 1;
            *SA
                .offset(
                    *fresh169 as isize,
                ) = if c0 as std::ffi::c_int <= c1 as std::ffi::c_int { p0 } else { t };
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        index = if p1 == 0 as std::ffi::c_int as std::ffi::c_long {
            i - 1 as std::ffi::c_int as std::ffi::c_long
        } else {
            index
        };
        *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let mut c0_0: uint16_t = *T
                .offset(
                    (p1
                        - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1_0: uint16_t = *T.offset(p1 as isize);
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c1_0 as sa_sint_t;
            let mut t_0: sa_sint_t = c0_0 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh170 = &mut (*induction_bucket.offset(c1_0 as isize));
            *fresh170 -= 1;
            *SA
                .offset(
                    *fresh170 as isize,
                ) = if c0_0 as std::ffi::c_int <= c1_0 as std::ffi::c_int {
                p1
            } else {
                t_0
            };
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        index = if p == 0 as std::ffi::c_int as std::ffi::c_long { i } else { index };
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let mut c0_1: uint16_t = *T
                .offset(
                    (p
                        - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1_1: uint16_t = *T.offset(p as isize);
            *SA.offset(i as isize) = c1_1 as sa_sint_t;
            let mut t_1: sa_sint_t = c0_1 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh171 = &mut (*induction_bucket.offset(c1_1 as isize));
            *fresh171 -= 1;
            *SA
                .offset(
                    *fresh171 as isize,
                ) = if c0_1 as std::ffi::c_int <= c1_1 as std::ffi::c_int {
                p
            } else {
                t_1
            };
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    index
}
unsafe extern "C" fn libsais16x64_final_bwt_aux_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut rm: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let mut c0: uint16_t = *T
                .offset(
                    (p0
                        - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1: uint16_t = *T.offset(p0 as isize);
            *SA
                .offset(
                    (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c1 as sa_sint_t;
            let mut t: sa_sint_t = c0 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh172 = &mut (*induction_bucket.offset(c1 as isize));
            *fresh172 -= 1;
            *SA
                .offset(
                    *fresh172 as isize,
                ) = if c0 as std::ffi::c_int <= c1 as std::ffi::c_int { p0 } else { t };
            if p0 & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p0 / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p0 as isize) as isize)
                    + 1 as std::ffi::c_int as std::ffi::c_long;
            }
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let mut c0_0: uint16_t = *T
                .offset(
                    (p1
                        - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1_0: uint16_t = *T.offset(p1 as isize);
            *SA
                .offset(
                    (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c1_0 as sa_sint_t;
            let mut t_0: sa_sint_t = c0_0 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh173 = &mut (*induction_bucket.offset(c1_0 as isize));
            *fresh173 -= 1;
            *SA
                .offset(
                    *fresh173 as isize,
                ) = if c0_0 as std::ffi::c_int <= c1_0 as std::ffi::c_int {
                p1
            } else {
                t_0
            };
            if p1 & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p1 / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p1 as isize) as isize)
                    + 1 as std::ffi::c_int as std::ffi::c_long;
            }
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let mut c0_1: uint16_t = *T
                .offset(
                    (p
                        - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                            as std::ffi::c_int as std::ffi::c_long) as isize,
                );
            let mut c1_1: uint16_t = *T.offset(p as isize);
            *SA.offset(i as isize) = c1_1 as sa_sint_t;
            let mut t_1: sa_sint_t = c0_1 as std::ffi::c_long | (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long);
            let fresh174 = &mut (*induction_bucket.offset(c1_1 as isize));
            *fresh174 -= 1;
            *SA
                .offset(
                    *fresh174 as isize,
                ) = if c0_1 as std::ffi::c_int <= c1_1 as std::ffi::c_int {
                p
            } else {
                t_1
            };
            if p & rm == 0 as std::ffi::c_int as std::ffi::c_long {
                *I
                    .offset(
                        (p / (rm + 1 as std::ffi::c_int as std::ffi::c_long)) as isize,
                    ) = *induction_bucket.offset(*T.offset(p as isize) as isize)
                    + 1 as std::ffi::c_int as std::ffi::c_long;
            }
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let fresh175 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            *fresh175 -= 1;
            *SA
                .offset(
                    *fresh175 as isize,
                ) = p0
                | (((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p0 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let fresh176 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            *fresh176 -= 1;
            *SA
                .offset(
                    *fresh176 as isize,
                ) = p1
                | (((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p1 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let fresh177 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            *fresh177 -= 1;
            *SA
                .offset(
                    *fresh177 as isize,
                ) = p
                | (((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_gsa_scan_right_to_left_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts0: *const uint16_t = (&*T.offset(s0 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts0 = Ts0.offset(-1);
        libsais16x64_prefetchr(
            (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts0
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - prefetch_distance - 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            );
        let mut Ts1: *const uint16_t = (&*T.offset(s1 as isize) as *const uint16_t)
            .offset(-(1 as std::ffi::c_int as isize));
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        Ts1 = Ts1.offset(-1);
        libsais16x64_prefetchr(
            (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                Ts1
            } else {
                std::ptr::null::<uint16_t>()
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long
            && *T.offset((p0 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int > 0 as std::ffi::c_int
        {
            p0 -= 1;
            let fresh178 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            *fresh178 -= 1;
            *SA
                .offset(
                    *fresh178 as isize,
                ) = p0
                | (((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p0 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long
            && *T.offset((p1 - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int > 0 as std::ffi::c_int
        {
            p1 -= 1;
            let fresh179 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            *fresh179 -= 1;
            *SA
                .offset(
                    *fresh179 as isize,
                ) = p1
                | (((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p1 as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long
            && *T.offset((p - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int > 0 as std::ffi::c_int
        {
            p -= 1;
            let fresh180 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            *fresh180 -= 1;
            *SA
                .offset(
                    *fresh180 as isize,
                ) = p
                | (((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) as std::ffi::c_int > *T.offset(p as isize) as std::ffi::c_int)
                    as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_right_to_left_32s(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchw(
            &mut *SA
                .offset(
                    (i - 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut s0: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts0: *const sa_sint_t = &*T
            .offset(
                (if s0 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s0
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts0.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s1: sa_sint_t = *SA
            .offset(
                (i - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        let mut Ts1: *const sa_sint_t = &*T
            .offset(
                (if s1 > 0 as std::ffi::c_int as std::ffi::c_long {
                    s1
                } else {
                    1 as std::ffi::c_int as std::ffi::c_long
                }) as isize,
            ) as *const sa_sint_t;
        libsais16x64_prefetchr(
            Ts1.offset(-(1 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
        );
        let mut s2: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s2 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s2 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s2 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut s3: sa_sint_t = *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
                    - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            );
        if s3 > 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_prefetchw(
                &mut *induction_bucket
                    .offset(
                        *T
                            .offset(
                                (s3 - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                            ) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
            );
            libsais16x64_prefetchr(
                (&*T.offset(s3 as isize) as *const sa_sint_t)
                    .offset(-(2 as std::ffi::c_int as isize)) as *const std::ffi::c_void,
            );
        }
        let mut p0: sa_sint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long;
        if p0 > 0 as std::ffi::c_int as std::ffi::c_long {
            p0 -= 1;
            let fresh181 = &mut (*induction_bucket
                .offset(*T.offset(p0 as isize) as isize));
            *fresh181 -= 1;
            *SA
                .offset(
                    *fresh181 as isize,
                ) = p0
                | (((*T
                    .offset(
                        (p0
                            - (p0 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) > *T.offset(p0 as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        let mut p1: sa_sint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        *SA
            .offset(
                (i - 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long;
        if p1 > 0 as std::ffi::c_int as std::ffi::c_long {
            p1 -= 1;
            let fresh182 = &mut (*induction_bucket
                .offset(*T.offset(p1 as isize) as isize));
            *fresh182 -= 1;
            *SA
                .offset(
                    *fresh182 as isize,
                ) = p1
                | (((*T
                    .offset(
                        (p1
                            - (p1 > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) > *T.offset(p1 as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 2 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        -= 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 1 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_sint_t = *SA.offset(i as isize);
        *SA.offset(i as isize) = p & 9223372036854775807 as std::ffi::c_long;
        if p > 0 as std::ffi::c_int as std::ffi::c_long {
            p -= 1;
            let fresh183 = &mut (*induction_bucket
                .offset(*T.offset(p as isize) as isize));
            *fresh183 -= 1;
            *SA
                .offset(
                    *fresh183 as isize,
                ) = p
                | (((*T
                    .offset(
                        (p
                            - (p > 0 as std::ffi::c_int as std::ffi::c_long)
                                as std::ffi::c_int as std::ffi::c_long) as isize,
                    ) > *T.offset(p as isize)) as std::ffi::c_int as sa_uint_t) << (64 as std::ffi::c_int - 1 as std::ffi::c_int)) as sa_sint_t;
        }
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_final_bwt_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut index: sa_sint_t = -(1 as std::ffi::c_int) as sa_sint_t;
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        index = libsais16x64_final_bwt_scan_right_to_left_16u(
            T,
            SA,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
    index
}
unsafe extern "C" fn libsais16x64_final_bwt_aux_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _k: sa_sint_t,
    mut rm: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_bwt_aux_scan_right_to_left_16u(
            T,
            SA,
            rm,
            I,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
    mut _k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || omp_block_size < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_sorting_scan_right_to_left_16u(
            T,
            SA,
            induction_bucket,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_gsa_scan_right_to_left_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
    mut _k: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || omp_block_size < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_gsa_scan_right_to_left_16u(
            T,
            SA,
            induction_bucket,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_final_sorting_scan_right_to_left_32s_omp(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut induction_bucket: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if threads == 1 as std::ffi::c_int as std::ffi::c_long
        || n < 65536 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_final_sorting_scan_right_to_left_32s(
            T,
            SA,
            induction_bucket,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
    }
}
unsafe extern "C" fn libsais16x64_clear_lms_suffixes_omp(
    mut SA: *mut sa_sint_t,
    mut _n: sa_sint_t,
    mut k: sa_sint_t,
    mut bucket_start: *mut sa_sint_t,
    mut bucket_end: *mut sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut c: fast_sint_t = 0;
    c = 0 as std::ffi::c_int as fast_sint_t;
    while c < k {
        if *bucket_end.offset(c as isize) > *bucket_start.offset(c as isize) {
            memset(
                &mut *SA.offset(*bucket_start.offset(c as isize) as isize)
                    as *mut sa_sint_t as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                (*bucket_end.offset(c as isize) as size_t)
                    .wrapping_sub(*bucket_start.offset(c as isize) as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        c += 1;
    }
}
unsafe extern "C" fn libsais16x64_induce_final_order_16u_omp(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut flags: sa_sint_t,
    mut r: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    if flags & 1 as std::ffi::c_int as std::ffi::c_long
        == 0 as std::ffi::c_int as std::ffi::c_long
    {
        if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
            *buckets
                .offset(
                    (6 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ) = *buckets
                .offset(
                    (7 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ) - 1 as std::ffi::c_int as std::ffi::c_long;
        }
        libsais16x64_final_sorting_scan_left_to_right_16u_omp(
            T,
            SA,
            n,
            k,
            &mut *buckets
                .offset(
                    (6 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            threads,
            thread_state,
        );
        if threads > 1 as std::ffi::c_int as std::ffi::c_long
            && n >= 65536 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_clear_lms_suffixes_omp(
                SA,
                n,
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as sa_sint_t,
                &mut *buckets
                    .offset(
                        (6 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                &mut *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
            );
        }
        if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0 {
            libsais16x64_flip_suffix_markers_omp(
                SA,
                *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
            );
            libsais16x64_final_gsa_scan_right_to_left_16u_omp(
                T,
                SA,
                *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                n
                    - *buckets
                        .offset(
                            (7 as std::ffi::c_int
                                * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                    << 8 as std::ffi::c_int)) as isize,
                        ),
                k,
                &mut *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
                thread_state,
            );
        } else {
            libsais16x64_final_sorting_scan_right_to_left_16u_omp(
                T,
                SA,
                0 as std::ffi::c_int as fast_sint_t,
                n,
                k,
                &mut *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
                thread_state,
            );
        }
        0 as std::ffi::c_int as sa_sint_t
    } else if !I.is_null() {
        libsais16x64_final_bwt_aux_scan_left_to_right_16u_omp(
            T,
            SA,
            n,
            k,
            r - 1 as std::ffi::c_int as std::ffi::c_long,
            I,
            &mut *buckets
                .offset(
                    (6 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            threads,
            thread_state,
        );
        if threads > 1 as std::ffi::c_int as std::ffi::c_long
            && n >= 65536 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_clear_lms_suffixes_omp(
                SA,
                n,
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as sa_sint_t,
                &mut *buckets
                    .offset(
                        (6 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                &mut *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
            );
        }
        libsais16x64_final_bwt_aux_scan_right_to_left_16u_omp(
            T,
            SA,
            n,
            k,
            r - 1 as std::ffi::c_int as std::ffi::c_long,
            I,
            &mut *buckets
                .offset(
                    (7 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            threads,
            thread_state,
        );
        return 0 as std::ffi::c_int as sa_sint_t;
    } else {
        libsais16x64_final_bwt_scan_left_to_right_16u_omp(
            T,
            SA,
            n,
            k,
            &mut *buckets
                .offset(
                    (6 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            threads,
            thread_state,
        );
        if threads > 1 as std::ffi::c_int as std::ffi::c_long
            && n >= 65536 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_clear_lms_suffixes_omp(
                SA,
                n,
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as sa_sint_t,
                &mut *buckets
                    .offset(
                        (6 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                &mut *buckets
                    .offset(
                        (7 as std::ffi::c_int
                            * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                                << 8 as std::ffi::c_int)) as isize,
                    ),
                threads,
            );
        }
        return libsais16x64_final_bwt_scan_right_to_left_16u_omp(
            T,
            SA,
            n,
            k,
            &mut *buckets
                .offset(
                    (7 as std::ffi::c_int
                        * (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                            << 8 as std::ffi::c_int)) as isize,
                ),
            threads,
            thread_state,
        );
    }
}
unsafe extern "C" fn libsais16x64_induce_final_order_32s_6k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((5 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_final_order_32s_4k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((2 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((3 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_final_order_32s_2k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((1 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(
        T,
        SA,
        n,
        &mut *buckets.offset((0 as std::ffi::c_int as std::ffi::c_long * k) as isize),
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_induce_final_order_32s_1k(
    mut T: *const sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_start_32s_1k(k, buckets);
    libsais16x64_final_sorting_scan_left_to_right_32s_omp(
        T,
        SA,
        n,
        buckets,
        threads,
        thread_state,
    );
    libsais16x64_count_suffixes_32s(T, n, k, buckets);
    libsais16x64_initialize_buckets_end_32s_1k(k, buckets);
    libsais16x64_final_sorting_scan_right_to_left_32s_omp(
        T,
        SA,
        n,
        buckets,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut f: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) -> sa_sint_t {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAm: *mut sa_sint_t = &mut *SA.offset(m as isize) as *mut sa_sint_t;
    let mut i: sa_sint_t = 0;
    let mut j: sa_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size
        - 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *SAm
                .offset(
                    (*SA
                        .offset(
                            (i
                                + 2 as std::ffi::c_int as std::ffi::c_long
                                    * prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as sa_uint_t >> 1 as std::ffi::c_int) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        let mut q0: sa_uint_t = *SA
            .offset(
                (i + prefetch_distance + 0 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) as sa_uint_t;
        let mut Tq0: *mut sa_sint_t = &mut *T.offset(q0 as isize) as *mut sa_sint_t;
        libsais16x64_prefetchw(
            (if *SAm.offset((q0 >> 1 as std::ffi::c_int) as isize)
                < 0 as std::ffi::c_int as std::ffi::c_long
            {
                Tq0
            } else {
                &mut *SAm.offset((q0 >> 1 as std::ffi::c_int) as isize) as *mut sa_sint_t
            }) as *const std::ffi::c_void,
        );
        let mut q1: sa_uint_t = *SA
            .offset(
                (i + prefetch_distance + 1 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) as sa_uint_t;
        let mut Tq1: *mut sa_sint_t = &mut *T.offset(q1 as isize) as *mut sa_sint_t;
        libsais16x64_prefetchw(
            (if *SAm.offset((q1 >> 1 as std::ffi::c_int) as isize)
                < 0 as std::ffi::c_int as std::ffi::c_long
            {
                Tq1
            } else {
                &mut *SAm.offset((q1 >> 1 as std::ffi::c_int) as isize) as *mut sa_sint_t
            }) as *const std::ffi::c_void,
        );
        let mut q2: sa_uint_t = *SA
            .offset(
                (i + prefetch_distance + 2 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) as sa_uint_t;
        let mut Tq2: *mut sa_sint_t = &mut *T.offset(q2 as isize) as *mut sa_sint_t;
        libsais16x64_prefetchw(
            (if *SAm.offset((q2 >> 1 as std::ffi::c_int) as isize)
                < 0 as std::ffi::c_int as std::ffi::c_long
            {
                Tq2
            } else {
                &mut *SAm.offset((q2 >> 1 as std::ffi::c_int) as isize) as *mut sa_sint_t
            }) as *const std::ffi::c_void,
        );
        let mut q3: sa_uint_t = *SA
            .offset(
                (i + prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long)
                    as isize,
            ) as sa_uint_t;
        let mut Tq3: *mut sa_sint_t = &mut *T.offset(q3 as isize) as *mut sa_sint_t;
        libsais16x64_prefetchw(
            (if *SAm.offset((q3 >> 1 as std::ffi::c_int) as isize)
                < 0 as std::ffi::c_int as std::ffi::c_long
            {
                Tq3
            } else {
                &mut *SAm.offset((q3 >> 1 as std::ffi::c_int) as isize) as *mut sa_sint_t
            }) as *const std::ffi::c_void,
        );
        let mut p0: sa_uint_t = *SA
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        let mut s0: sa_sint_t = *SAm.offset((p0 >> 1 as std::ffi::c_int) as isize);
        if s0 < 0 as std::ffi::c_int as std::ffi::c_long {
            let fresh184 = &mut (*T.offset(p0 as isize));
            *fresh184
                |= -(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long;
            f += 1;
            s0 = i + 0 as std::ffi::c_int as std::ffi::c_long
                + (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) + f;
        }
        *SAm.offset((p0 >> 1 as std::ffi::c_int) as isize) = s0 - f;
        let mut p1: sa_uint_t = *SA
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        let mut s1: sa_sint_t = *SAm.offset((p1 >> 1 as std::ffi::c_int) as isize);
        if s1 < 0 as std::ffi::c_int as std::ffi::c_long {
            let fresh185 = &mut (*T.offset(p1 as isize));
            *fresh185
                |= -(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long;
            f += 1;
            s1 = i + 1 as std::ffi::c_int as std::ffi::c_long
                + (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) + f;
        }
        *SAm.offset((p1 >> 1 as std::ffi::c_int) as isize) = s1 - f;
        let mut p2: sa_uint_t = *SA
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        let mut s2: sa_sint_t = *SAm.offset((p2 >> 1 as std::ffi::c_int) as isize);
        if s2 < 0 as std::ffi::c_int as std::ffi::c_long {
            let fresh186 = &mut (*T.offset(p2 as isize));
            *fresh186
                |= -(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long;
            f += 1;
            s2 = i + 2 as std::ffi::c_int as std::ffi::c_long
                + (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) + f;
        }
        *SAm.offset((p2 >> 1 as std::ffi::c_int) as isize) = s2 - f;
        let mut p3: sa_uint_t = *SA
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        let mut s3: sa_sint_t = *SAm.offset((p3 >> 1 as std::ffi::c_int) as isize);
        if s3 < 0 as std::ffi::c_int as std::ffi::c_long {
            let fresh187 = &mut (*T.offset(p3 as isize));
            *fresh187
                |= -(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long;
            f += 1;
            s3 = i + 3 as std::ffi::c_int as std::ffi::c_long
                + (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) + f;
        }
        *SAm.offset((p3 >> 1 as std::ffi::c_int) as isize) = s3 - f;
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j
        += 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance
            + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut p: sa_uint_t = *SA.offset(i as isize) as sa_uint_t;
        let mut s: sa_sint_t = *SAm.offset((p >> 1 as std::ffi::c_int) as isize);
        if s < 0 as std::ffi::c_int as std::ffi::c_long {
            let fresh188 = &mut (*T.offset(p as isize));
            *fresh188
                |= -(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long;
            f += 1;
            s = i
                + (-(9223372036854775807 as std::ffi::c_long)
                    - 1 as std::ffi::c_int as std::ffi::c_long) + f;
        }
        *SAm.offset((p >> 1 as std::ffi::c_int) as isize) = s - f;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    f
}
unsafe extern "C" fn libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s(
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut pl: *mut fast_sint_t,
    mut pr: *mut fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAl: *mut sa_uint_t = &mut *SA.offset(0 as std::ffi::c_int as isize)
        as *mut sa_sint_t as *mut sa_uint_t;
    let mut SAr: *mut sa_uint_t = &mut *SA.offset(0 as std::ffi::c_int as isize)
        as *mut sa_sint_t as *mut sa_uint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut l: fast_sint_t = *pl - 1 as std::ffi::c_int as std::ffi::c_long;
    let mut r: fast_sint_t = *pr - 1 as std::ffi::c_int as std::ffi::c_long;
    i = m + omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = m + omp_block_start + 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        libsais16x64_prefetchr(
            &mut *SA.offset((i - prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut p0: sa_uint_t = *SA
            .offset((i - 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SAl
            .offset(
                l as isize,
            ) = p0 & 9223372036854775807 as std::ffi::c_long as std::ffi::c_ulong;
        l
            -= ((p0 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        *SAr
            .offset(
                r as isize,
            ) = p0.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong);
        r
            -= (p0 as sa_sint_t > 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut p1: sa_uint_t = *SA
            .offset((i - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SAl
            .offset(
                l as isize,
            ) = p1 & 9223372036854775807 as std::ffi::c_long as std::ffi::c_ulong;
        l
            -= ((p1 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        *SAr
            .offset(
                r as isize,
            ) = p1.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong);
        r
            -= (p1 as sa_sint_t > 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut p2: sa_uint_t = *SA
            .offset((i - 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SAl
            .offset(
                l as isize,
            ) = p2 & 9223372036854775807 as std::ffi::c_long as std::ffi::c_ulong;
        l
            -= ((p2 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        *SAr
            .offset(
                r as isize,
            ) = p2.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong);
        r
            -= (p2 as sa_sint_t > 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        let mut p3: sa_uint_t = *SA
            .offset((i - 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            as sa_uint_t;
        *SAl
            .offset(
                l as isize,
            ) = p3 & 9223372036854775807 as std::ffi::c_long as std::ffi::c_ulong;
        l
            -= ((p3 as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        *SAr
            .offset(
                r as isize,
            ) = p3.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong);
        r
            -= (p3 as sa_sint_t > 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        i -= 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j -= 3 as std::ffi::c_int as std::ffi::c_long;
    while i >= j {
        let mut p: sa_uint_t = *SA.offset(i as isize) as sa_uint_t;
        *SAl
            .offset(
                l as isize,
            ) = p & 9223372036854775807 as std::ffi::c_long as std::ffi::c_ulong;
        l
            -= ((p as sa_sint_t) < 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        *SAr
            .offset(
                r as isize,
            ) = p.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong);
        r
            -= (p as sa_sint_t > 0 as std::ffi::c_int as std::ffi::c_long)
                as std::ffi::c_int as std::ffi::c_long;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
    *pl = l + 1 as std::ffi::c_int as std::ffi::c_long;
    *pr = r + 1 as std::ffi::c_int as std::ffi::c_long;
}
unsafe extern "C" fn libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut f: sa_sint_t = 0 as std::ffi::c_int as sa_sint_t;
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (m / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        m - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        f = libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s(
            T,
            SA,
            m,
            0 as std::ffi::c_int as sa_sint_t,
            omp_block_start,
            omp_block_size,
        );
    }
    f
}
unsafe extern "C" fn libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut f: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = ((n >> 1 as std::ffi::c_int) / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        (n >> 1 as std::ffi::c_int) - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        let mut l: fast_sint_t = m;
        let mut r: fast_sint_t = n + fs;
        libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s(
            SA,
            m,
            &mut l,
            &mut r,
            omp_block_start,
            omp_block_size,
        );
    }
    memcpy(
        &mut *SA.offset((n + fs - m) as isize) as *mut sa_sint_t
            as *mut std::ffi::c_void,
        &mut *SA.offset((m - f) as isize) as *mut sa_sint_t as *const std::ffi::c_void,
        (f as size_t)
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
    );
}
unsafe extern "C" fn libsais16x64_compact_lms_suffixes_32s_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut f: sa_sint_t = libsais16x64_renumber_unique_and_nonunique_lms_suffixes_32s_omp(
        T,
        SA,
        m,
        threads,
        thread_state,
    );
    libsais16x64_compact_unique_and_nonunique_lms_suffixes_32s_omp(
        SA,
        n,
        m,
        fs,
        f,
        threads,
        thread_state,
    );
    f
}
unsafe extern "C" fn libsais16x64_merge_unique_lms_suffixes_32s(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut l: fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAnm: *const sa_sint_t = &mut *SA
        .offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long + l) as isize)
        as *mut sa_sint_t;
    let mut i: sa_sint_t = 0;
    let mut j: sa_sint_t = 0;
    let fresh189 = SAnm;
    SAnm = SAnm.offset(1);
    let mut tmp: fast_sint_t = *fresh189;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - 6 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *T.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        let mut c0: sa_sint_t = *T
            .offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        if c0 < 0 as std::ffi::c_int as std::ffi::c_long {
            *T
                .offset(
                    (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c0 & 9223372036854775807 as std::ffi::c_long;
            *SA.offset(tmp as isize) = i + 0 as std::ffi::c_int as std::ffi::c_long;
            i += 1;
            let fresh190 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh190;
        }
        let mut c1: sa_sint_t = *T
            .offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        if c1 < 0 as std::ffi::c_int as std::ffi::c_long {
            *T
                .offset(
                    (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c1 & 9223372036854775807 as std::ffi::c_long;
            *SA.offset(tmp as isize) = i + 1 as std::ffi::c_int as std::ffi::c_long;
            i += 1;
            let fresh191 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh191;
        }
        let mut c2: sa_sint_t = *T
            .offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        if c2 < 0 as std::ffi::c_int as std::ffi::c_long {
            *T
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c2 & 9223372036854775807 as std::ffi::c_long;
            *SA.offset(tmp as isize) = i + 2 as std::ffi::c_int as std::ffi::c_long;
            i += 1;
            let fresh192 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh192;
        }
        let mut c3: sa_sint_t = *T
            .offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        if c3 < 0 as std::ffi::c_int as std::ffi::c_long {
            *T
                .offset(
                    (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                ) = c3 & 9223372036854775807 as std::ffi::c_long;
            *SA.offset(tmp as isize) = i + 3 as std::ffi::c_int as std::ffi::c_long;
            i += 1;
            let fresh193 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh193;
        }
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 6 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        let mut c: sa_sint_t = *T.offset(i as isize);
        if c < 0 as std::ffi::c_int as std::ffi::c_long {
            *T.offset(i as isize) = c & 9223372036854775807 as std::ffi::c_long;
            *SA.offset(tmp as isize) = i;
            i += 1;
            let fresh194 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh194;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_merge_nonunique_lms_suffixes_32s(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut l: fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut SAnm: *const sa_sint_t = &mut *SA
        .offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long + l) as isize)
        as *mut sa_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let fresh195 = SAnm;
    SAnm = SAnm.offset(1);
    let mut tmp: sa_sint_t = *fresh195;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *SA.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        if *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            == 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize) = tmp;
            let fresh196 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh196;
        }
        if *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            == 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize) = tmp;
            let fresh197 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh197;
        }
        if *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            == 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize) = tmp;
            let fresh198 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh198;
        }
        if *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            == 0 as std::ffi::c_int as std::ffi::c_long
        {
            *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize) = tmp;
            let fresh199 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh199;
        }
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        if *SA.offset(i as isize) == 0 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(i as isize) = tmp;
            let fresh200 = SAnm;
            SAnm = SAnm.offset(1);
            tmp = *fresh200;
        }
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_merge_unique_lms_suffixes_32s_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_merge_unique_lms_suffixes_32s(
            T,
            SA,
            n,
            m,
            0 as std::ffi::c_int as fast_sint_t,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_merge_nonunique_lms_suffixes_32s_omp(
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut f: sa_sint_t,
    mut _threads: sa_sint_t,
    mut _thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (m / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        m - omp_block_start
    };
    if omp_num_threads == 1 as std::ffi::c_int as std::ffi::c_long {
        libsais16x64_merge_nonunique_lms_suffixes_32s(
            SA,
            n,
            m,
            f,
            omp_block_start,
            omp_block_size,
        );
    }
}
unsafe extern "C" fn libsais16x64_merge_compacted_lms_suffixes_32s_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut f: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    libsais16x64_merge_unique_lms_suffixes_32s_omp(T, SA, n, m, threads, thread_state);
    libsais16x64_merge_nonunique_lms_suffixes_32s_omp(
        SA,
        n,
        m,
        f,
        threads,
        thread_state,
    );
}
unsafe extern "C" fn libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut f: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if f > 0 as std::ffi::c_int as std::ffi::c_long {
        memmove(
            &mut *SA.offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as *mut sa_sint_t as *mut std::ffi::c_void,
            &mut *SA.offset((n + fs - m) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
            (f as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_count_and_gather_compacted_lms_suffixes_32s_2k_omp(
            T,
            SA,
            n,
            k,
            buckets,
            threads,
            thread_state,
        );
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m - f, threads);
        memcpy(
            &mut *SA
                .offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long + f) as isize)
                as *mut sa_sint_t as *mut std::ffi::c_void,
            &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
            (m as size_t)
                .wrapping_sub(f as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        memset(
            &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
                as *mut std::ffi::c_void,
            0 as std::ffi::c_int,
            (m as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_merge_compacted_lms_suffixes_32s_omp(
            T,
            SA,
            n,
            m,
            f,
            threads,
            thread_state,
        );
    } else {
        libsais16x64_count_and_gather_lms_suffixes_32s_2k(
            T,
            SA,
            n,
            k,
            buckets,
            0 as std::ffi::c_int as fast_sint_t,
            n,
        );
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m, threads);
    };
}
unsafe extern "C" fn libsais16x64_reconstruct_compacted_lms_suffixes_32s_1k_omp(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut m: sa_sint_t,
    mut fs: sa_sint_t,
    mut f: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) {
    if f > 0 as std::ffi::c_int as std::ffi::c_long {
        memmove(
            &mut *SA.offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as *mut sa_sint_t as *mut std::ffi::c_void,
            &mut *SA.offset((n + fs - m) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
            (f as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_gather_compacted_lms_suffixes_32s(T, SA, n);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m - f, threads);
        memcpy(
            &mut *SA
                .offset((n - m - 1 as std::ffi::c_int as std::ffi::c_long + f) as isize)
                as *mut sa_sint_t as *mut std::ffi::c_void,
            &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
            (m as size_t)
                .wrapping_sub(f as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        memset(
            &mut *SA.offset(0 as std::ffi::c_int as isize) as *mut sa_sint_t
                as *mut std::ffi::c_void,
            0 as std::ffi::c_int,
            (m as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_merge_compacted_lms_suffixes_32s_omp(
            T,
            SA,
            n,
            m,
            f,
            threads,
            thread_state,
        );
    } else {
        libsais16x64_gather_lms_suffixes_32s(T, SA, n);
        libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m, threads);
    };
}
unsafe extern "C" fn libsais16x64_convert_32u_to_64u(
    mut S: *mut uint32_t,
    mut D: *mut uint64_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size;
    while i < j {
        *D.offset(i as isize) = *S.offset(i as isize) as uint64_t;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_convert_inplace_32u_to_64u(
    mut V: *mut uint32_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start + omp_block_size - 1 as std::ffi::c_int as std::ffi::c_long;
    j = omp_block_start;
    while i >= j {
        *V
            .offset(
                (i + i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *V.offset(i as isize);
        *V
            .offset(
                (i + i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = 0 as std::ffi::c_int as uint32_t;
        i -= 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_convert_inplace_64u_to_32u(
    mut V: *mut uint32_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size;
    while i < j {
        *V
            .offset(
                i as isize,
            ) = *V.offset((i + i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_convert_inplace_32u_to_64u_omp(
    mut V: *mut uint32_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    while n >= 65536 as std::ffi::c_int as std::ffi::c_long {
        let mut block_size: fast_sint_t = n >> 1 as std::ffi::c_int;
        n -= block_size;
        let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
        let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
        let mut omp_block_stride: fast_sint_t = (block_size / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
        let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
        let mut omp_block_size: fast_sint_t = if omp_thread_num
            < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
        {
            omp_block_stride
        } else {
            block_size - omp_block_start
        };
        libsais16x64_convert_32u_to_64u(
            (V as *mut std::ffi::c_void as *mut uint32_t).offset(n as isize),
            (V as *mut std::ffi::c_void as *mut uint64_t).offset(n as isize),
            omp_block_start,
            omp_block_size,
        );
    }
    libsais16x64_convert_inplace_32u_to_64u(V, 0 as std::ffi::c_int as fast_sint_t, n);
}
unsafe extern "C" fn libsais16x64_main_32s_recursion(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut fs: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
    mut local_buffer: *mut sa_sint_t,
) -> sa_sint_t {
    fs = if fs < 9223372036854775807 as std::ffi::c_long - n {
        fs
    } else {
        9223372036854775807 as std::ffi::c_long - n
    };
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long {
        let mut new_fs: sa_sint_t = if fs + fs + n + n
            <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        {
            fs + fs + n
        } else {
            2147483647 as std::ffi::c_int as std::ffi::c_long - n
        };
        if new_fs / k >= 6 as std::ffi::c_int as std::ffi::c_long
            || new_fs / k >= 4 as std::ffi::c_int as std::ffi::c_long
                && n
                    <= (2147483647 as std::ffi::c_int / 2 as std::ffi::c_int)
                        as std::ffi::c_long
            || new_fs / k < 4 as std::ffi::c_int as std::ffi::c_long && new_fs >= fs
        {
            libsais16x64_convert_inplace_64u_to_32u(
                T as *mut std::ffi::c_void as *mut uint32_t,
                0 as std::ffi::c_int as fast_sint_t,
                n,
            );
            let mut index: sa_sint_t = libsais16_int(
                T as *mut int32_t,
                SA as *mut int32_t,
                n as int32_t,
                k as int32_t,
                new_fs as int32_t,
            ) as sa_sint_t;
            if index >= 0 as std::ffi::c_int as std::ffi::c_long {
                libsais16x64_convert_inplace_32u_to_64u_omp(
                    SA as *mut uint32_t,
                    n,
                    threads,
                );
                libsais16x64_convert_inplace_32u_to_64u_omp(
                    T as *mut uint32_t,
                    n,
                    threads,
                );
            }
            return index;
        }
    }
    if k > 0 as std::ffi::c_int as std::ffi::c_long
        && (fs / k >= 6 as std::ffi::c_int as std::ffi::c_long
            || 1024 as std::ffi::c_int as std::ffi::c_long / k
                >= 6 as std::ffi::c_int as std::ffi::c_long
                && threads == 1 as std::ffi::c_int as std::ffi::c_long)
    {
        let mut alignment: sa_sint_t = if (fs
            - 1024 as std::ffi::c_int as std::ffi::c_long) / k
            >= 6 as std::ffi::c_int as std::ffi::c_long
        {
            1024 as std::ffi::c_int as sa_sint_t
        } else {
            16 as std::ffi::c_int as sa_sint_t
        };
        let mut buckets: *mut sa_sint_t = if (fs - alignment) / k
            >= 6 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_align_up(
                &mut *SA
                    .offset(
                        (n + fs - 6 as std::ffi::c_int as std::ffi::c_long * k
                            - alignment) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
                (alignment as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            ) as *mut sa_sint_t
        } else {
            &mut *SA
                .offset((n + fs - 6 as std::ffi::c_int as std::ffi::c_long * k) as isize)
                as *mut sa_sint_t
        };
        buckets = if 1024 as std::ffi::c_int as std::ffi::c_long / k
            >= 6 as std::ffi::c_int as std::ffi::c_long
            && threads == 1 as std::ffi::c_int as std::ffi::c_long
        {
            local_buffer
        } else {
            buckets
        };
        let mut m: sa_sint_t = libsais16x64_count_and_gather_lms_suffixes_32s_4k_omp(
            T,
            SA,
            n,
            k,
            buckets,
            threads,
            thread_state,
        );
        if m > 1 as std::ffi::c_int as std::ffi::c_long {
            memset(
                SA as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                (n as size_t)
                    .wrapping_sub(m as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
            let mut first_lms_suffix: sa_sint_t = *SA.offset((n - m) as isize);
            let mut left_suffixes_count: sa_sint_t = libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_6k(
                T,
                k,
                buckets,
                first_lms_suffix,
            );
            libsais16x64_radix_sort_lms_suffixes_32s_6k_omp(
                T,
                SA,
                n,
                m,
                &mut *buckets
                    .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize),
                threads,
                thread_state,
            );
            if (n / 8192 as std::ffi::c_int as std::ffi::c_long) < k {
                libsais16x64_radix_sort_set_markers_32s_6k_omp(
                    SA,
                    k,
                    &mut *buckets
                        .offset((4 as std::ffi::c_int as std::ffi::c_long * k) as isize),
                    threads,
                );
            }
            if threads > 1 as std::ffi::c_int as std::ffi::c_long
                && n >= 65536 as std::ffi::c_int as std::ffi::c_long
            {
                memset(
                    &mut *SA.offset((n - m) as isize) as *mut sa_sint_t
                        as *mut std::ffi::c_void,
                    0 as std::ffi::c_int,
                    (m as size_t)
                        .wrapping_mul(
                            ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                        ),
                );
            }
            libsais16x64_initialize_buckets_for_partial_sorting_32s_6k(
                T,
                k,
                buckets,
                first_lms_suffix,
                left_suffixes_count,
            );
            libsais16x64_induce_partial_order_32s_6k_omp(
                T,
                SA,
                n,
                k,
                buckets,
                first_lms_suffix,
                left_suffixes_count,
                threads,
                thread_state,
            );
            let mut names: sa_sint_t = if (n
                / 8192 as std::ffi::c_int as std::ffi::c_long) < k
            {
                libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(
                    SA,
                    n,
                    m,
                    threads,
                    thread_state,
                )
            } else {
                libsais16x64_renumber_and_gather_lms_suffixes_omp(
                    SA,
                    n,
                    m,
                    fs,
                    threads,
                    thread_state,
                )
            };
            if names < m {
                let mut f: sa_sint_t = if (n
                    / 8192 as std::ffi::c_int as std::ffi::c_long) < k
                {
                    libsais16x64_compact_lms_suffixes_32s_omp(
                        T,
                        SA,
                        n,
                        m,
                        fs,
                        threads,
                        thread_state,
                    )
                } else {
                    0 as std::ffi::c_int as std::ffi::c_long
                };
                if libsais16x64_main_32s_recursion(
                    SA
                        .offset(n as isize)
                        .offset(fs as isize)
                        .offset(-(m as isize))
                        .offset(f as isize),
                    SA,
                    m - f,
                    names - f,
                    fs + n - 2 as std::ffi::c_int as std::ffi::c_long * m + f,
                    threads,
                    thread_state,
                    local_buffer,
                ) != 0 as std::ffi::c_int as std::ffi::c_long
                {
                    return -(2 as std::ffi::c_int) as sa_sint_t;
                }
                libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(
                    T,
                    SA,
                    n,
                    k,
                    m,
                    fs,
                    f,
                    buckets,
                    threads,
                    thread_state,
                );
            } else {
                libsais16x64_count_lms_suffixes_32s_2k(T, n, k, buckets);
            }
            libsais16x64_initialize_buckets_start_and_end_32s_4k(k, buckets);
            libsais16x64_place_lms_suffixes_histogram_32s_4k(SA, n, k, m, buckets);
            libsais16x64_induce_final_order_32s_4k(
                T,
                SA,
                n,
                k,
                buckets,
                threads,
                thread_state,
            );
        } else {
            *SA
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *SA.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
            libsais16x64_initialize_buckets_start_and_end_32s_6k(k, buckets);
            libsais16x64_place_lms_suffixes_histogram_32s_6k(SA, n, k, m, buckets);
            libsais16x64_induce_final_order_32s_6k(
                T,
                SA,
                n,
                k,
                buckets,
                threads,
                thread_state,
            );
        }
        0 as std::ffi::c_int as sa_sint_t
    } else if k > 0 as std::ffi::c_int as std::ffi::c_long
        && n
            <= 9223372036854775807 as std::ffi::c_long
                / 2 as std::ffi::c_int as std::ffi::c_long
        && (fs / k >= 4 as std::ffi::c_int as std::ffi::c_long
            || 1024 as std::ffi::c_int as std::ffi::c_long / k
                >= 4 as std::ffi::c_int as std::ffi::c_long
                && threads == 1 as std::ffi::c_int as std::ffi::c_long)
    {
        let mut alignment_0: sa_sint_t = if (fs
            - 1024 as std::ffi::c_int as std::ffi::c_long) / k
            >= 4 as std::ffi::c_int as std::ffi::c_long
        {
            1024 as std::ffi::c_int as sa_sint_t
        } else {
            16 as std::ffi::c_int as sa_sint_t
        };
        let mut buckets_0: *mut sa_sint_t = if (fs - alignment_0) / k
            >= 4 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_align_up(
                &mut *SA
                    .offset(
                        (n + fs - 4 as std::ffi::c_int as std::ffi::c_long * k
                            - alignment_0) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
                (alignment_0 as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            ) as *mut sa_sint_t
        } else {
            &mut *SA
                .offset((n + fs - 4 as std::ffi::c_int as std::ffi::c_long * k) as isize)
                as *mut sa_sint_t
        };
        buckets_0 = if 1024 as std::ffi::c_int as std::ffi::c_long / k
            >= 4 as std::ffi::c_int as std::ffi::c_long
            && threads == 1 as std::ffi::c_int as std::ffi::c_long
        {
            local_buffer
        } else {
            buckets_0
        };
        let mut m_0: sa_sint_t = libsais16x64_count_and_gather_lms_suffixes_32s_2k_omp(
            T,
            SA,
            n,
            k,
            buckets_0,
            threads,
            thread_state,
        );
        if m_0 > 1 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_initialize_buckets_for_radix_and_partial_sorting_32s_4k(
                T,
                k,
                buckets_0,
                *SA.offset((n - m_0) as isize),
            );
            libsais16x64_radix_sort_lms_suffixes_32s_2k_omp(
                T,
                SA,
                n,
                m_0,
                &mut *buckets_0.offset(1 as std::ffi::c_int as isize),
                threads,
                thread_state,
            );
            libsais16x64_radix_sort_set_markers_32s_4k_omp(
                SA,
                k,
                &mut *buckets_0.offset(1 as std::ffi::c_int as isize),
                threads,
            );
            libsais16x64_place_lms_suffixes_interval_32s_4k(
                SA,
                n,
                k,
                m_0 - 1 as std::ffi::c_int as std::ffi::c_long,
                buckets_0,
            );
            libsais16x64_induce_partial_order_32s_4k_omp(
                T,
                SA,
                n,
                k,
                buckets_0,
                threads,
                thread_state,
            );
            let mut names_0: sa_sint_t = libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_4k_omp(
                SA,
                n,
                m_0,
                threads,
                thread_state,
            );
            if names_0 < m_0 {
                let mut f_0: sa_sint_t = libsais16x64_compact_lms_suffixes_32s_omp(
                    T,
                    SA,
                    n,
                    m_0,
                    fs,
                    threads,
                    thread_state,
                );
                if libsais16x64_main_32s_recursion(
                    SA
                        .offset(n as isize)
                        .offset(fs as isize)
                        .offset(-(m_0 as isize))
                        .offset(f_0 as isize),
                    SA,
                    m_0 - f_0,
                    names_0 - f_0,
                    fs + n - 2 as std::ffi::c_int as std::ffi::c_long * m_0 + f_0,
                    threads,
                    thread_state,
                    local_buffer,
                ) != 0 as std::ffi::c_int as std::ffi::c_long
                {
                    return -(2 as std::ffi::c_int) as sa_sint_t;
                }
                libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(
                    T,
                    SA,
                    n,
                    k,
                    m_0,
                    fs,
                    f_0,
                    buckets_0,
                    threads,
                    thread_state,
                );
            } else {
                libsais16x64_count_lms_suffixes_32s_2k(T, n, k, buckets_0);
            }
        } else {
            *SA
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *SA.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        }
        libsais16x64_initialize_buckets_start_and_end_32s_4k(k, buckets_0);
        libsais16x64_place_lms_suffixes_histogram_32s_4k(SA, n, k, m_0, buckets_0);
        libsais16x64_induce_final_order_32s_4k(
            T,
            SA,
            n,
            k,
            buckets_0,
            threads,
            thread_state,
        );
        return 0 as std::ffi::c_int as sa_sint_t;
    } else if k > 0 as std::ffi::c_int as std::ffi::c_long
        && (fs / k >= 2 as std::ffi::c_int as std::ffi::c_long
            || 1024 as std::ffi::c_int as std::ffi::c_long / k
                >= 2 as std::ffi::c_int as std::ffi::c_long
                && threads == 1 as std::ffi::c_int as std::ffi::c_long)
    {
        let mut alignment_1: sa_sint_t = if (fs
            - 1024 as std::ffi::c_int as std::ffi::c_long) / k
            >= 2 as std::ffi::c_int as std::ffi::c_long
        {
            1024 as std::ffi::c_int as sa_sint_t
        } else {
            16 as std::ffi::c_int as sa_sint_t
        };
        let mut buckets_1: *mut sa_sint_t = if (fs - alignment_1) / k
            >= 2 as std::ffi::c_int as std::ffi::c_long
        {
            libsais16x64_align_up(
                &mut *SA
                    .offset(
                        (n + fs - 2 as std::ffi::c_int as std::ffi::c_long * k
                            - alignment_1) as isize,
                    ) as *mut sa_sint_t as *const std::ffi::c_void,
                (alignment_1 as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            ) as *mut sa_sint_t
        } else {
            &mut *SA
                .offset((n + fs - 2 as std::ffi::c_int as std::ffi::c_long * k) as isize)
                as *mut sa_sint_t
        };
        buckets_1 = if 1024 as std::ffi::c_int as std::ffi::c_long / k
            >= 2 as std::ffi::c_int as std::ffi::c_long
            && threads == 1 as std::ffi::c_int as std::ffi::c_long
        {
            local_buffer
        } else {
            buckets_1
        };
        let mut m_1: sa_sint_t = libsais16x64_count_and_gather_lms_suffixes_32s_2k_omp(
            T,
            SA,
            n,
            k,
            buckets_1,
            threads,
            thread_state,
        );
        if m_1 > 1 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_32s_2k(
                T,
                k,
                buckets_1,
                *SA.offset((n - m_1) as isize),
            );
            libsais16x64_radix_sort_lms_suffixes_32s_2k_omp(
                T,
                SA,
                n,
                m_1,
                &mut *buckets_1.offset(1 as std::ffi::c_int as isize),
                threads,
                thread_state,
            );
            libsais16x64_place_lms_suffixes_interval_32s_2k(
                SA,
                n,
                k,
                m_1 - 1 as std::ffi::c_int as std::ffi::c_long,
                buckets_1,
            );
            libsais16x64_initialize_buckets_start_and_end_32s_2k(k, buckets_1);
            libsais16x64_induce_partial_order_32s_2k_omp(
                T,
                SA,
                n,
                k,
                buckets_1,
                threads,
                thread_state,
            );
            let mut names_1: sa_sint_t = libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(
                T,
                SA,
                n,
                m_1,
                threads,
            );
            if names_1 < m_1 {
                let mut f_1: sa_sint_t = libsais16x64_compact_lms_suffixes_32s_omp(
                    T,
                    SA,
                    n,
                    m_1,
                    fs,
                    threads,
                    thread_state,
                );
                if libsais16x64_main_32s_recursion(
                    SA
                        .offset(n as isize)
                        .offset(fs as isize)
                        .offset(-(m_1 as isize))
                        .offset(f_1 as isize),
                    SA,
                    m_1 - f_1,
                    names_1 - f_1,
                    fs + n - 2 as std::ffi::c_int as std::ffi::c_long * m_1 + f_1,
                    threads,
                    thread_state,
                    local_buffer,
                ) != 0 as std::ffi::c_int as std::ffi::c_long
                {
                    return -(2 as std::ffi::c_int) as sa_sint_t;
                }
                libsais16x64_reconstruct_compacted_lms_suffixes_32s_2k_omp(
                    T,
                    SA,
                    n,
                    k,
                    m_1,
                    fs,
                    f_1,
                    buckets_1,
                    threads,
                    thread_state,
                );
            } else {
                libsais16x64_count_lms_suffixes_32s_2k(T, n, k, buckets_1);
            }
        } else {
            *SA
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *SA.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        }
        libsais16x64_initialize_buckets_end_32s_2k(k, buckets_1);
        libsais16x64_place_lms_suffixes_histogram_32s_2k(SA, n, k, m_1, buckets_1);
        libsais16x64_initialize_buckets_start_and_end_32s_2k(k, buckets_1);
        libsais16x64_induce_final_order_32s_2k(
            T,
            SA,
            n,
            k,
            buckets_1,
            threads,
            thread_state,
        );
        return 0 as std::ffi::c_int as sa_sint_t;
    } else {
        let mut buffer: *mut sa_sint_t = if fs < k {
            libsais16x64_alloc_aligned(
                (k as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
                4096 as std::ffi::c_int as size_t,
            ) as *mut sa_sint_t
        } else {
            std::ptr::null_mut::<std::ffi::c_void>() as *mut sa_sint_t
        };
        let mut alignment_2: sa_sint_t = if fs
            - 1024 as std::ffi::c_int as std::ffi::c_long >= k
        {
            1024 as std::ffi::c_int as sa_sint_t
        } else {
            16 as std::ffi::c_int as sa_sint_t
        };
        let mut buckets_2: *mut sa_sint_t = if fs - alignment_2 >= k {
            libsais16x64_align_up(
                &mut *SA.offset((n + fs - k - alignment_2) as isize) as *mut sa_sint_t
                    as *const std::ffi::c_void,
                (alignment_2 as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            ) as *mut sa_sint_t
        } else if fs >= k {
            &mut *SA.offset((n + fs - k) as isize) as *mut sa_sint_t
        } else {
            buffer
        };
        if buckets_2.is_null() {
            return -(2 as std::ffi::c_int) as sa_sint_t;
        }
        memset(
            SA as *mut std::ffi::c_void,
            0 as std::ffi::c_int,
            (n as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_count_suffixes_32s(T, n, k, buckets_2);
        libsais16x64_initialize_buckets_end_32s_1k(k, buckets_2);
        let mut m_2: sa_sint_t = libsais16x64_radix_sort_lms_suffixes_32s_1k(
            T,
            SA,
            n,
            buckets_2,
        );
        if m_2 > 1 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_induce_partial_order_32s_1k_omp(
                T,
                SA,
                n,
                k,
                buckets_2,
                threads,
                thread_state,
            );
            let mut names_2: sa_sint_t = libsais16x64_renumber_and_mark_distinct_lms_suffixes_32s_1k_omp(
                T,
                SA,
                n,
                m_2,
                threads,
            );
            if names_2 < m_2 {
                if !buffer.is_null() {
                    libsais16x64_free_aligned(buffer as *mut std::ffi::c_void);
                    buckets_2 = std::ptr::null_mut::<sa_sint_t>();
                }
                let mut f_2: sa_sint_t = libsais16x64_compact_lms_suffixes_32s_omp(
                    T,
                    SA,
                    n,
                    m_2,
                    fs,
                    threads,
                    thread_state,
                );
                if libsais16x64_main_32s_recursion(
                    SA
                        .offset(n as isize)
                        .offset(fs as isize)
                        .offset(-(m_2 as isize))
                        .offset(f_2 as isize),
                    SA,
                    m_2 - f_2,
                    names_2 - f_2,
                    fs + n - 2 as std::ffi::c_int as std::ffi::c_long * m_2 + f_2,
                    threads,
                    thread_state,
                    local_buffer,
                ) != 0 as std::ffi::c_int as std::ffi::c_long
                {
                    return -(2 as std::ffi::c_int) as sa_sint_t;
                }
                libsais16x64_reconstruct_compacted_lms_suffixes_32s_1k_omp(
                    T,
                    SA,
                    n,
                    m_2,
                    fs,
                    f_2,
                    threads,
                    thread_state,
                );
                if buckets_2.is_null() {
                    buffer = libsais16x64_alloc_aligned(
                        (k as size_t)
                            .wrapping_mul(
                                ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                            ),
                        4096 as std::ffi::c_int as size_t,
                    ) as *mut sa_sint_t;
                    buckets_2 = buffer;
                }
                if buckets_2.is_null() {
                    return -(2 as std::ffi::c_int) as sa_sint_t;
                }
            }
            libsais16x64_count_suffixes_32s(T, n, k, buckets_2);
            libsais16x64_initialize_buckets_end_32s_1k(k, buckets_2);
            libsais16x64_place_lms_suffixes_interval_32s_1k(T, SA, k, m_2, buckets_2);
        }
        libsais16x64_induce_final_order_32s_1k(
            T,
            SA,
            n,
            k,
            buckets_2,
            threads,
            thread_state,
        );
        libsais16x64_free_aligned(buffer as *mut std::ffi::c_void);
        return 0 as std::ffi::c_int as sa_sint_t;
    }
}
unsafe extern "C" fn libsais16x64_main_32s_entry(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut fs: sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    let mut local_buffer: [sa_sint_t; 1024] = [0; 1024];
    libsais16x64_main_32s_recursion(
        T,
        SA,
        n,
        k,
        fs,
        threads,
        thread_state,
        local_buffer.as_mut_ptr(),
    )
}
unsafe extern "C" fn libsais16x64_main_16u(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut buckets: *mut sa_sint_t,
    mut flags: sa_sint_t,
    mut r: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut fs: sa_sint_t,
    mut freq: *mut sa_sint_t,
    mut threads: sa_sint_t,
    mut thread_state: *mut LIBSAIS_THREAD_STATE,
) -> sa_sint_t {
    fs = if fs < 9223372036854775807 as std::ffi::c_long - n {
        fs
    } else {
        9223372036854775807 as std::ffi::c_long - n
    };
    let mut m: sa_sint_t = libsais16x64_count_and_gather_lms_suffixes_16u_omp(
        T,
        SA,
        n,
        buckets,
        threads,
        thread_state,
    );
    let mut k: sa_sint_t = libsais16x64_initialize_buckets_start_and_end_16u(
        buckets,
        freq,
    );
    if flags & 2 as std::ffi::c_int as std::ffi::c_long != 0
        && (*buckets.offset(0 as std::ffi::c_int as isize)
            != 0 as std::ffi::c_int as std::ffi::c_long
            || *buckets.offset(2 as std::ffi::c_int as isize)
                != 0 as std::ffi::c_int as std::ffi::c_long
            || *buckets.offset(3 as std::ffi::c_int as isize)
                != 1 as std::ffi::c_int as std::ffi::c_long)
    {
        return -(1 as std::ffi::c_int) as sa_sint_t;
    }
    if m > 0 as std::ffi::c_int as std::ffi::c_long {
        let mut first_lms_suffix: sa_sint_t = *SA.offset((n - m) as isize);
        let mut left_suffixes_count: sa_sint_t = libsais16x64_initialize_buckets_for_lms_suffixes_radix_sort_16u(
            T,
            buckets,
            first_lms_suffix,
        );
        if threads > 1 as std::ffi::c_int as std::ffi::c_long
            && n >= 65536 as std::ffi::c_int as std::ffi::c_long
        {
            memset(
                SA as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                (n as size_t)
                    .wrapping_sub(m as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        libsais16x64_radix_sort_lms_suffixes_16u_omp(
            T,
            SA,
            n,
            m,
            flags,
            buckets,
            threads,
            thread_state,
        );
        if threads > 1 as std::ffi::c_int as std::ffi::c_long
            && n >= 65536 as std::ffi::c_int as std::ffi::c_long
        {
            memset(
                &mut *SA.offset((n - m) as isize) as *mut sa_sint_t
                    as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                (m as size_t)
                    .wrapping_mul(
                        ::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong,
                    ),
            );
        }
        libsais16x64_initialize_buckets_for_partial_sorting_16u(
            T,
            buckets,
            first_lms_suffix,
            left_suffixes_count,
        );
        libsais16x64_induce_partial_order_16u_omp(
            T,
            SA,
            n,
            k,
            flags,
            buckets,
            first_lms_suffix,
            left_suffixes_count,
            threads,
            thread_state,
        );
        let mut names: sa_sint_t = libsais16x64_renumber_and_gather_lms_suffixes_omp(
            SA,
            n,
            m,
            fs,
            threads,
            thread_state,
        );
        if names < m {
            if libsais16x64_main_32s_entry(
                SA.offset(n as isize).offset(fs as isize).offset(-(m as isize)),
                SA,
                m,
                names,
                fs + n - 2 as std::ffi::c_int as std::ffi::c_long * m,
                threads,
                thread_state,
            ) != 0 as std::ffi::c_int as std::ffi::c_long
            {
                return -(2 as std::ffi::c_int) as sa_sint_t;
            }
            libsais16x64_gather_lms_suffixes_16u_omp(T, SA, n, threads, thread_state);
            libsais16x64_reconstruct_lms_suffixes_omp(SA, n, m, threads);
        }
        libsais16x64_place_lms_suffixes_interval_16u(SA, n, m, flags, buckets);
    } else {
        memset(
            SA as *mut std::ffi::c_void,
            0 as std::ffi::c_int,
            (n as size_t)
                .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        );
    }
    libsais16x64_induce_final_order_16u_omp(
        T,
        SA,
        n,
        k,
        flags,
        r,
        I,
        buckets,
        threads,
        thread_state,
    )
}
unsafe extern "C" fn libsais16x64_main(
    mut T: *const uint16_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut flags: sa_sint_t,
    mut r: sa_sint_t,
    mut I: *mut sa_sint_t,
    mut fs: sa_sint_t,
    mut freq: *mut sa_sint_t,
    mut threads: sa_sint_t,
) -> sa_sint_t {
    let mut thread_state: *mut LIBSAIS_THREAD_STATE = if threads
        > 1 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_alloc_thread_state(threads)
    } else {
        std::ptr::null_mut::<LIBSAIS_THREAD_STATE>()
    };
    let mut buckets: *mut sa_sint_t = libsais16x64_alloc_aligned(
        (8 as std::ffi::c_int as size_t)
            .wrapping_mul(
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong,
            )
            .wrapping_mul(::core::mem::size_of::<sa_sint_t>() as std::ffi::c_ulong),
        4096 as std::ffi::c_int as size_t,
    ) as *mut sa_sint_t;
    let mut index: sa_sint_t = if !buckets.is_null()
        && (!thread_state.is_null()
            || threads == 1 as std::ffi::c_int as std::ffi::c_long)
    {
        libsais16x64_main_16u(
            T,
            SA,
            n,
            buckets,
            flags,
            r,
            I,
            fs,
            freq,
            threads,
            thread_state,
        )
    } else {
        -(2 as std::ffi::c_int) as std::ffi::c_long
    };
    libsais16x64_free_aligned(buckets as *mut std::ffi::c_void);
    libsais16x64_free_thread_state(thread_state);
    index
}
unsafe extern "C" fn libsais16x64_main_long(
    mut T: *mut sa_sint_t,
    mut SA: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut k: sa_sint_t,
    mut fs: sa_sint_t,
    mut threads: sa_sint_t,
) -> sa_sint_t {
    let mut thread_state: *mut LIBSAIS_THREAD_STATE = if threads
        > 1 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_alloc_thread_state(threads)
    } else {
        std::ptr::null_mut::<LIBSAIS_THREAD_STATE>()
    };
    let mut index: sa_sint_t = if !thread_state.is_null()
        || threads == 1 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_main_32s_entry(T, SA, n, k, fs, threads, thread_state)
    } else {
        -(2 as std::ffi::c_int) as std::ffi::c_long
    };
    libsais16x64_free_thread_state(thread_state);
    index
}
unsafe extern "C" fn libsais16x64_bwt_copy_16u(
    mut U: *mut uint16_t,
    mut A: *mut sa_sint_t,
    mut n: sa_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = 0 as std::ffi::c_int as fast_sint_t;
    j = n - 7 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &mut *A.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        *U
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 4 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 4 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 5 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 5 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 6 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 6 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        *U
            .offset(
                (i + 7 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *A.offset((i + 7 as std::ffi::c_int as std::ffi::c_long) as isize)
            as uint16_t;
        i += 8 as std::ffi::c_int as std::ffi::c_long;
    }
    j += 7 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *U.offset(i as isize) = *A.offset(i as isize) as uint16_t;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_bwt_copy_16u_omp(
    mut U: *mut uint16_t,
    mut A: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_block_start: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_size: fast_sint_t = n;
    libsais16x64_bwt_copy_16u(
        U.offset(omp_block_start as isize),
        A.offset(omp_block_start as isize),
        omp_block_size,
    );
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64(
    mut T: *const uint16_t,
    mut SA: *mut int64_t,
    mut n: int64_t,
    mut fs: int64_t,
    mut freq: *mut int64_t,
) -> int64_t {
    if T.is_null() || SA.is_null() || n < 0 as std::ffi::c_int as std::ffi::c_long
        || fs < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if !freq.is_null() {
            memset(
                freq as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong)
                    .wrapping_mul(::core::mem::size_of::<int64_t>() as std::ffi::c_ulong),
            );
        }
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(0 as std::ffi::c_int as isize) = 0 as std::ffi::c_int as int64_t;
            if !freq.is_null() {
                let fresh201 = &mut (*freq
                    .offset(*T.offset(0 as std::ffi::c_int as isize) as isize));
                *fresh201 += 1;
            }
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long {
        let mut new_fs: sa_sint_t = if fs + fs + n + n
            <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        {
            fs + fs + n
        } else {
            2147483647 as std::ffi::c_int as std::ffi::c_long - n
        };
        let mut index: sa_sint_t = libsais16(
            T,
            SA as *mut int32_t,
            n as int32_t,
            new_fs as int32_t,
            freq as *mut int32_t,
        ) as sa_sint_t;
        if index >= 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_convert_inplace_32u_to_64u_omp(
                SA as *mut uint32_t,
                n,
                1 as std::ffi::c_int as sa_sint_t,
            );
            if !freq.is_null() {
                libsais16x64_convert_inplace_32u_to_64u_omp(
                    freq as *mut uint32_t,
                    (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int) as sa_sint_t,
                    1 as std::ffi::c_int as sa_sint_t,
                );
            }
        }
        return index;
    }
    libsais16x64_main(
        T,
        SA,
        n,
        0 as std::ffi::c_int as sa_sint_t,
        0 as std::ffi::c_int as sa_sint_t,
        std::ptr::null_mut::<sa_sint_t>(),
        fs,
        freq,
        1 as std::ffi::c_int as sa_sint_t,
    )
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_gsa(
    mut T: *const uint16_t,
    mut SA: *mut int64_t,
    mut n: int64_t,
    mut fs: int64_t,
    mut freq: *mut int64_t,
) -> int64_t {
    if T.is_null() || SA.is_null() || n < 0 as std::ffi::c_int as std::ffi::c_long
        || n > 0 as std::ffi::c_int as std::ffi::c_long
            && *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int != 0 as std::ffi::c_int
        || fs < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if !freq.is_null() {
            memset(
                freq as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong)
                    .wrapping_mul(::core::mem::size_of::<int64_t>() as std::ffi::c_ulong),
            );
        }
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(0 as std::ffi::c_int as isize) = 0 as std::ffi::c_int as int64_t;
            if !freq.is_null() {
                let fresh202 = &mut (*freq
                    .offset(*T.offset(0 as std::ffi::c_int as isize) as isize));
                *fresh202 += 1;
            }
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long {
        let mut new_fs: sa_sint_t = if fs + fs + n + n
            <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        {
            fs + fs + n
        } else {
            2147483647 as std::ffi::c_int as std::ffi::c_long - n
        };
        let mut index: sa_sint_t = libsais16_gsa(
            T,
            SA as *mut int32_t,
            n as int32_t,
            new_fs as int32_t,
            freq as *mut int32_t,
        ) as sa_sint_t;
        if index >= 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_convert_inplace_32u_to_64u_omp(
                SA as *mut uint32_t,
                n,
                1 as std::ffi::c_int as sa_sint_t,
            );
            if !freq.is_null() {
                libsais16x64_convert_inplace_32u_to_64u_omp(
                    freq as *mut uint32_t,
                    (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int) as sa_sint_t,
                    1 as std::ffi::c_int as sa_sint_t,
                );
            }
        }
        return index;
    }
    libsais16x64_main(
        T,
        SA,
        n,
        2 as std::ffi::c_int as sa_sint_t,
        0 as std::ffi::c_int as sa_sint_t,
        std::ptr::null_mut::<sa_sint_t>(),
        fs,
        freq,
        1 as std::ffi::c_int as sa_sint_t,
    )
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_long(
    mut T: *mut int64_t,
    mut SA: *mut int64_t,
    mut n: int64_t,
    mut k: int64_t,
    mut fs: int64_t,
) -> int64_t {
    if T.is_null() || SA.is_null() || n < 0 as std::ffi::c_int as std::ffi::c_long
        || fs < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *SA.offset(0 as std::ffi::c_int as isize) = 0 as std::ffi::c_int as int64_t;
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    libsais16x64_main_long(T, SA, n, k, fs, 1 as std::ffi::c_int as sa_sint_t)
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_bwt(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut A: *mut int64_t,
    mut n: int64_t,
    mut fs: int64_t,
    mut freq: *mut int64_t,
) -> int64_t {
    if T.is_null() || U.is_null() || A.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
        || fs < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if !freq.is_null() {
            memset(
                freq as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong)
                    .wrapping_mul(::core::mem::size_of::<int64_t>() as std::ffi::c_ulong),
            );
        }
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *U
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *T.offset(0 as std::ffi::c_int as isize);
            if !freq.is_null() {
                let fresh203 = &mut (*freq
                    .offset(*T.offset(0 as std::ffi::c_int as isize) as isize));
                *fresh203 += 1;
            }
        }
        return n;
    }
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long {
        let mut new_fs: sa_sint_t = if fs + fs + n + n
            <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        {
            fs + fs + n
        } else {
            2147483647 as std::ffi::c_int as std::ffi::c_long - n
        };
        let mut index: sa_sint_t = libsais16_bwt(
            T,
            U,
            A as *mut int32_t,
            n as int32_t,
            new_fs as int32_t,
            freq as *mut int32_t,
        ) as sa_sint_t;
        if index >= 0 as std::ffi::c_int as std::ffi::c_long && !freq.is_null() {
            libsais16x64_convert_inplace_32u_to_64u_omp(
                freq as *mut uint32_t,
                (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as sa_sint_t,
                1 as std::ffi::c_int as sa_sint_t,
            );
        }
        return index;
    }
    let mut index_0: sa_sint_t = libsais16x64_main(
        T,
        A,
        n,
        1 as std::ffi::c_int as sa_sint_t,
        0 as std::ffi::c_int as sa_sint_t,
        std::ptr::null_mut::<sa_sint_t>(),
        fs,
        freq,
        1 as std::ffi::c_int as sa_sint_t,
    );
    if index_0 >= 0 as std::ffi::c_int as std::ffi::c_long {
        index_0 += 1;
        *U
            .offset(
                0 as std::ffi::c_int as isize,
            ) = *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        libsais16x64_bwt_copy_16u_omp(
            U.offset(1 as std::ffi::c_int as isize),
            A,
            index_0 - 1 as std::ffi::c_int as std::ffi::c_long,
            1 as std::ffi::c_int as sa_sint_t,
        );
        libsais16x64_bwt_copy_16u_omp(
            U.offset(index_0 as isize),
            A.offset(index_0 as isize),
            n - index_0,
            1 as std::ffi::c_int as sa_sint_t,
        );
    }
    index_0
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_bwt_aux(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut A: *mut int64_t,
    mut n: int64_t,
    mut fs: int64_t,
    mut freq: *mut int64_t,
    mut r: int64_t,
    mut I: *mut int64_t,
) -> int64_t {
    if T.is_null() || U.is_null() || A.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
        || fs < 0 as std::ffi::c_int as std::ffi::c_long
        || r < 2 as std::ffi::c_int as std::ffi::c_long
        || r & (r - 1 as std::ffi::c_int as std::ffi::c_long)
            != 0 as std::ffi::c_int as std::ffi::c_long || I.is_null()
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if !freq.is_null() {
            memset(
                freq as *mut std::ffi::c_void,
                0 as std::ffi::c_int,
                ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                    << 8 as std::ffi::c_int) as std::ffi::c_ulong)
                    .wrapping_mul(::core::mem::size_of::<int64_t>() as std::ffi::c_ulong),
            );
        }
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *U
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *T.offset(0 as std::ffi::c_int as isize);
            if !freq.is_null() {
                let fresh204 = &mut (*freq
                    .offset(*T.offset(0 as std::ffi::c_int as isize) as isize));
                *fresh204 += 1;
            }
        }
        *I.offset(0 as std::ffi::c_int as isize) = n;
        return 0 as std::ffi::c_int as int64_t;
    }
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        && r <= 2147483647 as std::ffi::c_int as std::ffi::c_long
    {
        let mut new_fs: sa_sint_t = if fs + fs + n + n
            <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        {
            fs + fs + n
        } else {
            2147483647 as std::ffi::c_int as std::ffi::c_long - n
        };
        let mut index: sa_sint_t = libsais16_bwt_aux(
            T,
            U,
            A as *mut int32_t,
            n as int32_t,
            new_fs as int32_t,
            freq as *mut int32_t,
            r as int32_t,
            I as *mut int32_t,
        ) as sa_sint_t;
        if index >= 0 as std::ffi::c_int as std::ffi::c_long {
            libsais16x64_convert_inplace_32u_to_64u_omp(
                I as *mut uint32_t,
                1 as std::ffi::c_int as std::ffi::c_long
                    + (n - 1 as std::ffi::c_int as std::ffi::c_long) / r,
                1 as std::ffi::c_int as sa_sint_t,
            );
            if !freq.is_null() {
                libsais16x64_convert_inplace_32u_to_64u_omp(
                    freq as *mut uint32_t,
                    (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int) as sa_sint_t,
                    1 as std::ffi::c_int as sa_sint_t,
                );
            }
        }
        return index;
    }
    let mut index_0: sa_sint_t = libsais16x64_main(
        T,
        A,
        n,
        1 as std::ffi::c_int as sa_sint_t,
        r,
        I,
        fs,
        freq,
        1 as std::ffi::c_int as sa_sint_t,
    );
    if index_0 == 0 as std::ffi::c_int as std::ffi::c_long {
        *U
            .offset(
                0 as std::ffi::c_int as isize,
            ) = *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        libsais16x64_bwt_copy_16u_omp(
            U.offset(1 as std::ffi::c_int as isize),
            A,
            *I.offset(0 as std::ffi::c_int as isize)
                - 1 as std::ffi::c_int as std::ffi::c_long,
            1 as std::ffi::c_int as sa_sint_t,
        );
        libsais16x64_bwt_copy_16u_omp(
            U.offset(*I.offset(0 as std::ffi::c_int as isize) as isize),
            A.offset(*I.offset(0 as std::ffi::c_int as isize) as isize),
            n - *I.offset(0 as std::ffi::c_int as isize),
            1 as std::ffi::c_int as sa_sint_t,
        );
    }
    index_0
}
unsafe extern "C" fn libsais16x64_unbwt_compute_histogram(
    mut T: *const uint16_t,
    mut n: fast_sint_t,
    mut count: *mut sa_uint_t,
) {
    let mut i: fast_sint_t = 0;
    i = 0 as std::ffi::c_int as fast_sint_t;
    while i < n {
        let fresh205 = &mut (*count.offset(*T.offset(i as isize) as isize));
        *fresh205 = (*fresh205).wrapping_add(1);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_unbwt_calculate_fastbits(
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
) {
    let mut v: fast_uint_t = 0;
    let mut w: fast_uint_t = 0;
    let mut sum: fast_uint_t = 0;
    v = 0 as std::ffi::c_int as fast_uint_t;
    sum = 1 as std::ffi::c_int as fast_uint_t;
    w = 0 as std::ffi::c_int as fast_uint_t;
    while w
        < (((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
            as std::ffi::c_ulong
    {
        let mut prev: fast_uint_t = sum;
        sum = (sum as std::ffi::c_ulong).wrapping_add(*bucket2.offset(w as isize))
            as fast_uint_t as fast_uint_t;
        *bucket2.offset(w as isize) = prev;
        if prev != sum {
            while v
                <= sum.wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong) >> shift
            {
                *fastbits.offset(v as isize) = w as uint16_t;
                v = v.wrapping_add(1);
            }
        }
        w = w.wrapping_add(1);
    }
}
unsafe extern "C" fn libsais16x64_unbwt_calculate_P(
    mut T: *const uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut index: fast_uint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_end: fast_sint_t,
) {
    let mut i: fast_sint_t = omp_block_start;
    let mut j: fast_sint_t = index as fast_sint_t;
    if omp_block_end < j {
        j = omp_block_end;
    }
    while i < j {
        let mut c: fast_uint_t = *T.offset(i as isize) as fast_uint_t;
        let fresh206 = &mut (*bucket2.offset(c as isize));
        let fresh207 = *fresh206;
        *fresh206 = (*fresh206).wrapping_add(1);
        *P.offset(fresh207 as isize) = i as sa_uint_t;
        i += 1;
    }
    let mut i_0: fast_sint_t = index as fast_sint_t;
    let mut j_0: fast_sint_t = omp_block_end;
    if omp_block_start > i_0 {
        i_0 = omp_block_start;
    }
    T = T.offset(-(1 as std::ffi::c_int as isize));
    i_0 += 1 as std::ffi::c_int as std::ffi::c_long;
    while i_0 <= j_0 {
        let mut c_0: fast_uint_t = *T.offset(i_0 as isize) as fast_uint_t;
        let fresh208 = &mut (*bucket2.offset(c_0 as isize));
        let fresh209 = *fresh208;
        *fresh208 = (*fresh208).wrapping_add(1);
        *P.offset(fresh209 as isize) = i_0 as sa_uint_t;
        i_0 += 1;
    }
}
unsafe extern "C" fn libsais16x64_unbwt_init_single(
    mut T: *const uint16_t,
    mut P: *mut sa_uint_t,
    mut n: sa_sint_t,
    mut freq: *const sa_sint_t,
    mut I: *const sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
) {
    let mut index: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
    let mut shift: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    while n >> shift > (1 as std::ffi::c_int as sa_sint_t) << 17 as std::ffi::c_int {
        shift = shift.wrapping_add(1);
    }
    if !freq.is_null() {
        memcpy(
            bucket2 as *mut std::ffi::c_void,
            freq as *const std::ffi::c_void,
            ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
                as std::ffi::c_ulong)
                .wrapping_mul(::core::mem::size_of::<sa_uint_t>() as std::ffi::c_ulong),
        );
    } else {
        memset(
            bucket2 as *mut std::ffi::c_void,
            0 as std::ffi::c_int,
            ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
                as std::ffi::c_ulong)
                .wrapping_mul(::core::mem::size_of::<sa_uint_t>() as std::ffi::c_ulong),
        );
        libsais16x64_unbwt_compute_histogram(T, n, bucket2);
    }
    libsais16x64_unbwt_calculate_fastbits(bucket2, fastbits, shift);
    libsais16x64_unbwt_calculate_P(
        T,
        P,
        bucket2,
        index,
        0 as std::ffi::c_int as fast_sint_t,
        n,
    );
}
unsafe extern "C" fn libsais16x64_unbwt_decode_1(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_2(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_3(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_4(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut i3: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut U3: *mut uint16_t = U2.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    let mut p3: fast_uint_t = *i3;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        let mut c3: uint16_t = *fastbits.offset((p3 >> shift) as isize);
        if *bucket2.offset(c3 as isize) <= p3 {
            loop {
                c3 = c3.wrapping_add(1);
                if *bucket2.offset(c3 as isize) > p3 {
                    break;
                }
            }
        }
        p3 = *P.offset(p3 as isize);
        *U3.offset(i as isize) = c3;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
    *i3 = p3;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_5(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut i3: *mut fast_uint_t,
    mut i4: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut U3: *mut uint16_t = U2.offset(r as isize);
    let mut U4: *mut uint16_t = U3.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    let mut p3: fast_uint_t = *i3;
    let mut p4: fast_uint_t = *i4;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        let mut c3: uint16_t = *fastbits.offset((p3 >> shift) as isize);
        if *bucket2.offset(c3 as isize) <= p3 {
            loop {
                c3 = c3.wrapping_add(1);
                if *bucket2.offset(c3 as isize) > p3 {
                    break;
                }
            }
        }
        p3 = *P.offset(p3 as isize);
        *U3.offset(i as isize) = c3;
        let mut c4: uint16_t = *fastbits.offset((p4 >> shift) as isize);
        if *bucket2.offset(c4 as isize) <= p4 {
            loop {
                c4 = c4.wrapping_add(1);
                if *bucket2.offset(c4 as isize) > p4 {
                    break;
                }
            }
        }
        p4 = *P.offset(p4 as isize);
        *U4.offset(i as isize) = c4;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
    *i3 = p3;
    *i4 = p4;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_6(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut i3: *mut fast_uint_t,
    mut i4: *mut fast_uint_t,
    mut i5: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut U3: *mut uint16_t = U2.offset(r as isize);
    let mut U4: *mut uint16_t = U3.offset(r as isize);
    let mut U5: *mut uint16_t = U4.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    let mut p3: fast_uint_t = *i3;
    let mut p4: fast_uint_t = *i4;
    let mut p5: fast_uint_t = *i5;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        let mut c3: uint16_t = *fastbits.offset((p3 >> shift) as isize);
        if *bucket2.offset(c3 as isize) <= p3 {
            loop {
                c3 = c3.wrapping_add(1);
                if *bucket2.offset(c3 as isize) > p3 {
                    break;
                }
            }
        }
        p3 = *P.offset(p3 as isize);
        *U3.offset(i as isize) = c3;
        let mut c4: uint16_t = *fastbits.offset((p4 >> shift) as isize);
        if *bucket2.offset(c4 as isize) <= p4 {
            loop {
                c4 = c4.wrapping_add(1);
                if *bucket2.offset(c4 as isize) > p4 {
                    break;
                }
            }
        }
        p4 = *P.offset(p4 as isize);
        *U4.offset(i as isize) = c4;
        let mut c5: uint16_t = *fastbits.offset((p5 >> shift) as isize);
        if *bucket2.offset(c5 as isize) <= p5 {
            loop {
                c5 = c5.wrapping_add(1);
                if *bucket2.offset(c5 as isize) > p5 {
                    break;
                }
            }
        }
        p5 = *P.offset(p5 as isize);
        *U5.offset(i as isize) = c5;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
    *i3 = p3;
    *i4 = p4;
    *i5 = p5;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_7(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut i3: *mut fast_uint_t,
    mut i4: *mut fast_uint_t,
    mut i5: *mut fast_uint_t,
    mut i6: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut U3: *mut uint16_t = U2.offset(r as isize);
    let mut U4: *mut uint16_t = U3.offset(r as isize);
    let mut U5: *mut uint16_t = U4.offset(r as isize);
    let mut U6: *mut uint16_t = U5.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    let mut p3: fast_uint_t = *i3;
    let mut p4: fast_uint_t = *i4;
    let mut p5: fast_uint_t = *i5;
    let mut p6: fast_uint_t = *i6;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        let mut c3: uint16_t = *fastbits.offset((p3 >> shift) as isize);
        if *bucket2.offset(c3 as isize) <= p3 {
            loop {
                c3 = c3.wrapping_add(1);
                if *bucket2.offset(c3 as isize) > p3 {
                    break;
                }
            }
        }
        p3 = *P.offset(p3 as isize);
        *U3.offset(i as isize) = c3;
        let mut c4: uint16_t = *fastbits.offset((p4 >> shift) as isize);
        if *bucket2.offset(c4 as isize) <= p4 {
            loop {
                c4 = c4.wrapping_add(1);
                if *bucket2.offset(c4 as isize) > p4 {
                    break;
                }
            }
        }
        p4 = *P.offset(p4 as isize);
        *U4.offset(i as isize) = c4;
        let mut c5: uint16_t = *fastbits.offset((p5 >> shift) as isize);
        if *bucket2.offset(c5 as isize) <= p5 {
            loop {
                c5 = c5.wrapping_add(1);
                if *bucket2.offset(c5 as isize) > p5 {
                    break;
                }
            }
        }
        p5 = *P.offset(p5 as isize);
        *U5.offset(i as isize) = c5;
        let mut c6: uint16_t = *fastbits.offset((p6 >> shift) as isize);
        if *bucket2.offset(c6 as isize) <= p6 {
            loop {
                c6 = c6.wrapping_add(1);
                if *bucket2.offset(c6 as isize) > p6 {
                    break;
                }
            }
        }
        p6 = *P.offset(p6 as isize);
        *U6.offset(i as isize) = c6;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
    *i3 = p3;
    *i4 = p4;
    *i5 = p5;
    *i6 = p6;
}
unsafe extern "C" fn libsais16x64_unbwt_decode_8(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut shift: fast_uint_t,
    mut r: fast_uint_t,
    mut i0: *mut fast_uint_t,
    mut i1: *mut fast_uint_t,
    mut i2: *mut fast_uint_t,
    mut i3: *mut fast_uint_t,
    mut i4: *mut fast_uint_t,
    mut i5: *mut fast_uint_t,
    mut i6: *mut fast_uint_t,
    mut i7: *mut fast_uint_t,
    mut k: fast_uint_t,
) {
    let mut U0: *mut uint16_t = U;
    let mut U1: *mut uint16_t = U0.offset(r as isize);
    let mut U2: *mut uint16_t = U1.offset(r as isize);
    let mut U3: *mut uint16_t = U2.offset(r as isize);
    let mut U4: *mut uint16_t = U3.offset(r as isize);
    let mut U5: *mut uint16_t = U4.offset(r as isize);
    let mut U6: *mut uint16_t = U5.offset(r as isize);
    let mut U7: *mut uint16_t = U6.offset(r as isize);
    let mut i: fast_uint_t = 0;
    let mut p0: fast_uint_t = *i0;
    let mut p1: fast_uint_t = *i1;
    let mut p2: fast_uint_t = *i2;
    let mut p3: fast_uint_t = *i3;
    let mut p4: fast_uint_t = *i4;
    let mut p5: fast_uint_t = *i5;
    let mut p6: fast_uint_t = *i6;
    let mut p7: fast_uint_t = *i7;
    i = 0 as std::ffi::c_int as fast_uint_t;
    while i != k {
        let mut c0: uint16_t = *fastbits.offset((p0 >> shift) as isize);
        if *bucket2.offset(c0 as isize) <= p0 {
            loop {
                c0 = c0.wrapping_add(1);
                if *bucket2.offset(c0 as isize) > p0 {
                    break;
                }
            }
        }
        p0 = *P.offset(p0 as isize);
        *U0.offset(i as isize) = c0;
        let mut c1: uint16_t = *fastbits.offset((p1 >> shift) as isize);
        if *bucket2.offset(c1 as isize) <= p1 {
            loop {
                c1 = c1.wrapping_add(1);
                if *bucket2.offset(c1 as isize) > p1 {
                    break;
                }
            }
        }
        p1 = *P.offset(p1 as isize);
        *U1.offset(i as isize) = c1;
        let mut c2: uint16_t = *fastbits.offset((p2 >> shift) as isize);
        if *bucket2.offset(c2 as isize) <= p2 {
            loop {
                c2 = c2.wrapping_add(1);
                if *bucket2.offset(c2 as isize) > p2 {
                    break;
                }
            }
        }
        p2 = *P.offset(p2 as isize);
        *U2.offset(i as isize) = c2;
        let mut c3: uint16_t = *fastbits.offset((p3 >> shift) as isize);
        if *bucket2.offset(c3 as isize) <= p3 {
            loop {
                c3 = c3.wrapping_add(1);
                if *bucket2.offset(c3 as isize) > p3 {
                    break;
                }
            }
        }
        p3 = *P.offset(p3 as isize);
        *U3.offset(i as isize) = c3;
        let mut c4: uint16_t = *fastbits.offset((p4 >> shift) as isize);
        if *bucket2.offset(c4 as isize) <= p4 {
            loop {
                c4 = c4.wrapping_add(1);
                if *bucket2.offset(c4 as isize) > p4 {
                    break;
                }
            }
        }
        p4 = *P.offset(p4 as isize);
        *U4.offset(i as isize) = c4;
        let mut c5: uint16_t = *fastbits.offset((p5 >> shift) as isize);
        if *bucket2.offset(c5 as isize) <= p5 {
            loop {
                c5 = c5.wrapping_add(1);
                if *bucket2.offset(c5 as isize) > p5 {
                    break;
                }
            }
        }
        p5 = *P.offset(p5 as isize);
        *U5.offset(i as isize) = c5;
        let mut c6: uint16_t = *fastbits.offset((p6 >> shift) as isize);
        if *bucket2.offset(c6 as isize) <= p6 {
            loop {
                c6 = c6.wrapping_add(1);
                if *bucket2.offset(c6 as isize) > p6 {
                    break;
                }
            }
        }
        p6 = *P.offset(p6 as isize);
        *U6.offset(i as isize) = c6;
        let mut c7: uint16_t = *fastbits.offset((p7 >> shift) as isize);
        if *bucket2.offset(c7 as isize) <= p7 {
            loop {
                c7 = c7.wrapping_add(1);
                if *bucket2.offset(c7 as isize) > p7 {
                    break;
                }
            }
        }
        p7 = *P.offset(p7 as isize);
        *U7.offset(i as isize) = c7;
        i = i.wrapping_add(1);
    }
    *i0 = p0;
    *i1 = p1;
    *i2 = p2;
    *i3 = p3;
    *i4 = p4;
    *i5 = p5;
    *i6 = p6;
    *i7 = p7;
}
unsafe extern "C" fn libsais16x64_unbwt_decode(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut n: sa_sint_t,
    mut r: sa_sint_t,
    mut I: *const sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut blocks: fast_sint_t,
    mut remainder: fast_uint_t,
) {
    let mut shift: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    while n >> shift > (1 as std::ffi::c_int as sa_sint_t) << 17 as std::ffi::c_int {
        shift = shift.wrapping_add(1);
    }
    let mut offset: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    while blocks > 8 as std::ffi::c_int as std::ffi::c_long {
        let mut i0: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        let mut i4: fast_uint_t = *I.offset(4 as std::ffi::c_int as isize);
        let mut i5: fast_uint_t = *I.offset(5 as std::ffi::c_int as isize);
        let mut i6: fast_uint_t = *I.offset(6 as std::ffi::c_int as isize);
        let mut i7: fast_uint_t = *I.offset(7 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_8(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0,
            &mut i1,
            &mut i2,
            &mut i3,
            &mut i4,
            &mut i5,
            &mut i6,
            &mut i7,
            r as fast_uint_t,
        );
        I = I.offset(8 as std::ffi::c_int as isize);
        blocks -= 8 as std::ffi::c_int as std::ffi::c_long;
        offset = (offset as std::ffi::c_ulong)
            .wrapping_add(
                (8 as std::ffi::c_int as std::ffi::c_ulong)
                    .wrapping_mul(r as fast_uint_t),
            ) as fast_uint_t as fast_uint_t;
    }
    if blocks == 1 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_0: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_1(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            &mut i0_0,
            remainder,
        );
    } else if blocks == 2 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_1: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_0: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_2(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_1,
            &mut i1_0,
            remainder,
        );
        libsais16x64_unbwt_decode_1(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            &mut i0_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else if blocks == 3 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_2: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_1: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_0: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_3(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_2,
            &mut i1_1,
            &mut i2_0,
            remainder,
        );
        libsais16x64_unbwt_decode_2(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_2,
            &mut i1_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else if blocks == 4 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_3: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_2: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_1: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3_0: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_4(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_3,
            &mut i1_2,
            &mut i2_1,
            &mut i3_0,
            remainder,
        );
        libsais16x64_unbwt_decode_3(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_3,
            &mut i1_2,
            &mut i2_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else if blocks == 5 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_4: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_3: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_2: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3_1: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        let mut i4_0: fast_uint_t = *I.offset(4 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_5(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_4,
            &mut i1_3,
            &mut i2_2,
            &mut i3_1,
            &mut i4_0,
            remainder,
        );
        libsais16x64_unbwt_decode_4(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_4,
            &mut i1_3,
            &mut i2_2,
            &mut i3_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else if blocks == 6 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_5: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_4: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_3: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3_2: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        let mut i4_1: fast_uint_t = *I.offset(4 as std::ffi::c_int as isize);
        let mut i5_0: fast_uint_t = *I.offset(5 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_6(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_5,
            &mut i1_4,
            &mut i2_3,
            &mut i3_2,
            &mut i4_1,
            &mut i5_0,
            remainder,
        );
        libsais16x64_unbwt_decode_5(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_5,
            &mut i1_4,
            &mut i2_3,
            &mut i3_2,
            &mut i4_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else if blocks == 7 as std::ffi::c_int as std::ffi::c_long {
        let mut i0_6: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_5: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_4: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3_3: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        let mut i4_2: fast_uint_t = *I.offset(4 as std::ffi::c_int as isize);
        let mut i5_1: fast_uint_t = *I.offset(5 as std::ffi::c_int as isize);
        let mut i6_0: fast_uint_t = *I.offset(6 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_7(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_6,
            &mut i1_5,
            &mut i2_4,
            &mut i3_3,
            &mut i4_2,
            &mut i5_1,
            &mut i6_0,
            remainder,
        );
        libsais16x64_unbwt_decode_6(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_6,
            &mut i1_5,
            &mut i2_4,
            &mut i3_3,
            &mut i4_2,
            &mut i5_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    } else {
        let mut i0_7: fast_uint_t = *I.offset(0 as std::ffi::c_int as isize);
        let mut i1_6: fast_uint_t = *I.offset(1 as std::ffi::c_int as isize);
        let mut i2_5: fast_uint_t = *I.offset(2 as std::ffi::c_int as isize);
        let mut i3_4: fast_uint_t = *I.offset(3 as std::ffi::c_int as isize);
        let mut i4_3: fast_uint_t = *I.offset(4 as std::ffi::c_int as isize);
        let mut i5_2: fast_uint_t = *I.offset(5 as std::ffi::c_int as isize);
        let mut i6_1: fast_uint_t = *I.offset(6 as std::ffi::c_int as isize);
        let mut i7_0: fast_uint_t = *I.offset(7 as std::ffi::c_int as isize);
        libsais16x64_unbwt_decode_8(
            U.offset(offset as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_7,
            &mut i1_6,
            &mut i2_5,
            &mut i3_4,
            &mut i4_3,
            &mut i5_2,
            &mut i6_1,
            &mut i7_0,
            remainder,
        );
        libsais16x64_unbwt_decode_7(
            U.offset(offset as isize).offset(remainder as isize),
            P,
            bucket2,
            fastbits,
            shift,
            r as fast_uint_t,
            &mut i0_7,
            &mut i1_6,
            &mut i2_5,
            &mut i3_4,
            &mut i4_3,
            &mut i5_2,
            &mut i6_1,
            (r as fast_uint_t).wrapping_sub(remainder),
        );
    };
}
unsafe extern "C" fn libsais16x64_unbwt_decode_omp(
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut n: sa_sint_t,
    mut r: sa_sint_t,
    mut I: *const sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut _threads: sa_sint_t,
) {
    let mut blocks: fast_sint_t = 1 as std::ffi::c_int as std::ffi::c_long
        + (n - 1 as std::ffi::c_int as std::ffi::c_long) / r;
    let mut remainder: fast_uint_t = (n as fast_uint_t)
        .wrapping_sub(
            (r as fast_uint_t)
                .wrapping_mul(
                    (blocks as fast_uint_t)
                        .wrapping_sub(1 as std::ffi::c_int as std::ffi::c_ulong),
                ),
        );
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = blocks / omp_num_threads;
    let mut omp_block_remainder: fast_sint_t = blocks % omp_num_threads;
    let mut omp_block_size: fast_sint_t = omp_block_stride
        + (omp_thread_num < omp_block_remainder) as std::ffi::c_int as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_block_stride * omp_thread_num
        + (if omp_thread_num < omp_block_remainder {
            omp_thread_num
        } else {
            omp_block_remainder
        });
    libsais16x64_unbwt_decode(
        U.offset((r * omp_block_start) as isize),
        P,
        n,
        r,
        I.offset(omp_block_start as isize),
        bucket2,
        fastbits,
        omp_block_size,
        if omp_thread_num < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long {
            r as fast_uint_t
        } else {
            remainder
        },
    );
}
unsafe extern "C" fn libsais16x64_unbwt_core(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut n: sa_sint_t,
    mut freq: *const sa_sint_t,
    mut r: sa_sint_t,
    mut I: *const sa_uint_t,
    mut bucket2: *mut sa_uint_t,
    mut fastbits: *mut uint16_t,
    mut _buckets: *mut sa_uint_t,
    mut threads: sa_sint_t,
) -> sa_sint_t {
    libsais16x64_unbwt_init_single(T, P, n, freq, I, bucket2, fastbits);
    libsais16x64_unbwt_decode_omp(U, P, n, r, I, bucket2, fastbits, threads);
    0 as std::ffi::c_int as sa_sint_t
}
unsafe extern "C" fn libsais16x64_unbwt_main(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut P: *mut sa_uint_t,
    mut n: sa_sint_t,
    mut freq: *const sa_sint_t,
    mut r: sa_sint_t,
    mut I: *const sa_uint_t,
    mut threads: sa_sint_t,
) -> sa_sint_t {
    let mut shift: fast_uint_t = 0 as std::ffi::c_int as fast_uint_t;
    while n >> shift > (1 as std::ffi::c_int as sa_sint_t) << 17 as std::ffi::c_int {
        shift = shift.wrapping_add(1);
    }
    let mut bucket2: *mut sa_uint_t = libsais16x64_alloc_aligned(
        ((((1 as std::ffi::c_int) << 8 as std::ffi::c_int) << 8 as std::ffi::c_int)
            as std::ffi::c_ulong)
            .wrapping_mul(::core::mem::size_of::<sa_uint_t>() as std::ffi::c_ulong),
        4096 as std::ffi::c_int as size_t,
    ) as *mut sa_uint_t;
    let mut fastbits: *mut uint16_t = libsais16x64_alloc_aligned(
        (1 as std::ffi::c_int as size_t)
            .wrapping_add((n >> shift) as size_t)
            .wrapping_mul(::core::mem::size_of::<uint16_t>() as std::ffi::c_ulong),
        4096 as std::ffi::c_int as size_t,
    ) as *mut uint16_t;
    let mut buckets: *mut sa_uint_t = if threads
        > 1 as std::ffi::c_int as std::ffi::c_long
        && n >= 262144 as std::ffi::c_int as std::ffi::c_long
    {
        libsais16x64_alloc_aligned(
            (threads as size_t)
                .wrapping_mul(
                    (((1 as std::ffi::c_int) << 8 as std::ffi::c_int)
                        << 8 as std::ffi::c_int) as std::ffi::c_ulong,
                )
                .wrapping_mul(::core::mem::size_of::<sa_uint_t>() as std::ffi::c_ulong),
            4096 as std::ffi::c_int as size_t,
        ) as *mut sa_uint_t
    } else {
        std::ptr::null_mut::<sa_uint_t>()
    };
    let mut index: sa_sint_t = if !bucket2.is_null() && !fastbits.is_null()
        && (!buckets.is_null() || threads == 1 as std::ffi::c_int as std::ffi::c_long
            || n < 262144 as std::ffi::c_int as std::ffi::c_long)
    {
        libsais16x64_unbwt_core(
            T,
            U,
            P,
            n,
            freq,
            r,
            I,
            bucket2,
            fastbits,
            buckets,
            threads,
        )
    } else {
        -(2 as std::ffi::c_int) as std::ffi::c_long
    };
    libsais16x64_free_aligned(buckets as *mut std::ffi::c_void);
    libsais16x64_free_aligned(fastbits as *mut std::ffi::c_void);
    libsais16x64_free_aligned(bucket2 as *mut std::ffi::c_void);
    index
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_unbwt(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut A: *mut int64_t,
    mut n: int64_t,
    mut freq: *const int64_t,
    mut i: int64_t,
) -> int64_t {
    libsais16x64_unbwt_aux(T, U, A, n, freq, n, &i)
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_unbwt_aux(
    mut T: *const uint16_t,
    mut U: *mut uint16_t,
    mut A: *mut int64_t,
    mut n: int64_t,
    mut freq: *const int64_t,
    mut r: int64_t,
    mut I: *const int64_t,
) -> int64_t {
    if T.is_null() || U.is_null() || A.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
        || r != n
            && (r < 2 as std::ffi::c_int as std::ffi::c_long
                || r & (r - 1 as std::ffi::c_int as std::ffi::c_long)
                    != 0 as std::ffi::c_int as std::ffi::c_long) || I.is_null()
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if *I.offset(0 as std::ffi::c_int as isize) != n {
            return -(1 as std::ffi::c_int) as int64_t;
        }
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *U
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *T.offset(0 as std::ffi::c_int as isize);
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    let mut t: fast_sint_t = 0;
    t = 0 as std::ffi::c_int as fast_sint_t;
    while t <= (n - 1 as std::ffi::c_int as std::ffi::c_long) / r {
        if *I.offset(t as isize) <= 0 as std::ffi::c_int as std::ffi::c_long
            || *I.offset(t as isize) > n
        {
            return -(1 as std::ffi::c_int) as int64_t;
        }
        t += 1;
    }
    if n <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        && r <= 2147483647 as std::ffi::c_int as std::ffi::c_long
        && (n - 1 as std::ffi::c_int as std::ffi::c_long) / r
            < 1024 as std::ffi::c_int as std::ffi::c_long
    {
        let mut indexes: [int32_t; 1024] = [0; 1024];
        t = 0 as std::ffi::c_int as fast_sint_t;
        while t <= (n - 1 as std::ffi::c_int as std::ffi::c_long) / r {
            indexes[t as usize] = *I.offset(t as isize) as int32_t;
            t += 1;
        }
        return libsais16_unbwt_aux(
            T,
            U,
            A as *mut int32_t,
            n as int32_t,
            std::ptr::null::<int32_t>(),
            r as int32_t,
            indexes.as_mut_ptr(),
        ) as int64_t;
    }
    libsais16x64_unbwt_main(
        T,
        U,
        A as *mut sa_uint_t,
        n,
        freq,
        r,
        I as *const sa_uint_t,
        1 as std::ffi::c_int as sa_sint_t,
    )
}
unsafe extern "C" fn libsais16x64_compute_phi(
    mut SA: *const sa_sint_t,
    mut PLCP: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut k: sa_sint_t = if omp_block_start > 0 as std::ffi::c_int as std::ffi::c_long
    {
        *SA.offset((omp_block_start - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
    } else {
        n
    };
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &*SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        *PLCP
            .offset(
                *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ) = k;
        k = *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize);
        *PLCP
            .offset(
                *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ) = k;
        k = *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize);
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        *PLCP
            .offset(
                *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ) = k;
        k = *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize);
        *PLCP
            .offset(
                *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            ) = k;
        k = *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize);
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *PLCP.offset(*SA.offset(i as isize) as isize) = k;
        k = *SA.offset(i as isize);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_compute_phi_omp(
    mut SA: *const sa_sint_t,
    mut PLCP: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    libsais16x64_compute_phi(SA, PLCP, n, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_compute_plcp(
    mut T: *const uint16_t,
    mut PLCP: *mut sa_sint_t,
    mut n: fast_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut l: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance;
    while i < j {
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T.offset((*PLCP.offset((i + prefetch_distance) as isize) + l) as isize)
                as *const uint16_t as *const std::ffi::c_void,
        );
        let mut k: fast_sint_t = *PLCP.offset(i as isize);
        let mut m: fast_sint_t = n - (if i > k { i } else { k });
        while l < m
            && *T.offset((i + l) as isize) as std::ffi::c_int
                == *T.offset((k + l) as isize) as std::ffi::c_int
        {
            l += 1;
        }
        *PLCP.offset(i as isize) = l;
        l
            -= (l != 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance;
    while i < j {
        let mut k_0: fast_sint_t = *PLCP.offset(i as isize);
        let mut m_0: fast_sint_t = n - (if i > k_0 { i } else { k_0 });
        while l < m_0
            && *T.offset((i + l) as isize) as std::ffi::c_int
                == *T.offset((k_0 + l) as isize) as std::ffi::c_int
        {
            l += 1;
        }
        *PLCP.offset(i as isize) = l;
        l
            -= (l != 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_compute_plcp_omp(
    mut T: *const uint16_t,
    mut PLCP: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    libsais16x64_compute_plcp(T, PLCP, n, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_compute_plcp_gsa(
    mut T: *const uint16_t,
    mut PLCP: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    let mut l: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance;
    while i < j {
        libsais16x64_prefetchw(
            &mut *PLCP
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *mut sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*T.offset((*PLCP.offset((i + prefetch_distance) as isize) + l) as isize)
                as *const uint16_t as *const std::ffi::c_void,
        );
        let mut k: fast_sint_t = *PLCP.offset(i as isize);
        while *T.offset((i + l) as isize) as std::ffi::c_int > 0 as std::ffi::c_int
            && *T.offset((i + l) as isize) as std::ffi::c_int
                == *T.offset((k + l) as isize) as std::ffi::c_int
        {
            l += 1;
        }
        *PLCP.offset(i as isize) = l;
        l
            -= (l != 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance;
    while i < j {
        let mut k_0: fast_sint_t = *PLCP.offset(i as isize);
        while *T.offset((i + l) as isize) as std::ffi::c_int > 0 as std::ffi::c_int
            && *T.offset((i + l) as isize) as std::ffi::c_int
                == *T.offset((k_0 + l) as isize) as std::ffi::c_int
        {
            l += 1;
        }
        *PLCP.offset(i as isize) = l;
        l
            -= (l != 0 as std::ffi::c_int as std::ffi::c_long) as std::ffi::c_int
                as std::ffi::c_long;
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_compute_plcp_gsa_omp(
    mut T: *const uint16_t,
    mut PLCP: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    libsais16x64_compute_plcp_gsa(T, PLCP, omp_block_start, omp_block_size);
}
unsafe extern "C" fn libsais16x64_compute_lcp(
    mut PLCP: *const sa_sint_t,
    mut SA: *const sa_sint_t,
    mut LCP: *mut sa_sint_t,
    mut omp_block_start: fast_sint_t,
    mut omp_block_size: fast_sint_t,
) {
    let prefetch_distance: fast_sint_t = 32 as std::ffi::c_int as fast_sint_t;
    let mut i: fast_sint_t = 0;
    let mut j: fast_sint_t = 0;
    i = omp_block_start;
    j = omp_block_start + omp_block_size - prefetch_distance
        - 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        libsais16x64_prefetchr(
            &*SA
                .offset(
                    (i + 2 as std::ffi::c_int as std::ffi::c_long * prefetch_distance)
                        as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchw(
            &mut *LCP.offset((i + prefetch_distance) as isize) as *mut sa_sint_t
                as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        *LCP
            .offset(
                (i + 0 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *PLCP
            .offset(
                *SA.offset((i + 0 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        *LCP
            .offset(
                (i + 1 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *PLCP
            .offset(
                *SA.offset((i + 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        libsais16x64_prefetchr(
            &*PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        libsais16x64_prefetchr(
            &*PLCP
                .offset(
                    *SA
                        .offset(
                            (i + prefetch_distance
                                + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
                        ) as isize,
                ) as *const sa_sint_t as *const std::ffi::c_void,
        );
        *LCP
            .offset(
                (i + 2 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *PLCP
            .offset(
                *SA.offset((i + 2 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        *LCP
            .offset(
                (i + 3 as std::ffi::c_int as std::ffi::c_long) as isize,
            ) = *PLCP
            .offset(
                *SA.offset((i + 3 as std::ffi::c_int as std::ffi::c_long) as isize)
                    as isize,
            );
        i += 4 as std::ffi::c_int as std::ffi::c_long;
    }
    j += prefetch_distance + 3 as std::ffi::c_int as std::ffi::c_long;
    while i < j {
        *LCP.offset(i as isize) = *PLCP.offset(*SA.offset(i as isize) as isize);
        i += 1 as std::ffi::c_int as std::ffi::c_long;
    }
}
unsafe extern "C" fn libsais16x64_compute_lcp_omp(
    mut PLCP: *const sa_sint_t,
    mut SA: *const sa_sint_t,
    mut LCP: *mut sa_sint_t,
    mut n: sa_sint_t,
    mut _threads: sa_sint_t,
) {
    let mut omp_thread_num: fast_sint_t = 0 as std::ffi::c_int as fast_sint_t;
    let mut omp_num_threads: fast_sint_t = 1 as std::ffi::c_int as fast_sint_t;
    let mut omp_block_stride: fast_sint_t = (n / omp_num_threads) & -(16 as std::ffi::c_int) as std::ffi::c_long;
    let mut omp_block_start: fast_sint_t = omp_thread_num * omp_block_stride;
    let mut omp_block_size: fast_sint_t = if omp_thread_num
        < omp_num_threads - 1 as std::ffi::c_int as std::ffi::c_long
    {
        omp_block_stride
    } else {
        n - omp_block_start
    };
    libsais16x64_compute_lcp(PLCP, SA, LCP, omp_block_start, omp_block_size);
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_plcp(
    mut T: *const uint16_t,
    mut SA: *const int64_t,
    mut PLCP: *mut int64_t,
    mut n: int64_t,
) -> int64_t {
    if T.is_null() || SA.is_null() || PLCP.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *PLCP
                .offset(0 as std::ffi::c_int as isize) = 0 as std::ffi::c_int as int64_t;
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    libsais16x64_compute_phi_omp(SA, PLCP, n, 1 as std::ffi::c_int as sa_sint_t);
    libsais16x64_compute_plcp_omp(T, PLCP, n, 1 as std::ffi::c_int as sa_sint_t);
    0 as std::ffi::c_int as int64_t
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_plcp_gsa(
    mut T: *const uint16_t,
    mut SA: *const int64_t,
    mut PLCP: *mut int64_t,
    mut n: int64_t,
) -> int64_t {
    if T.is_null() || SA.is_null() || PLCP.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
        || n > 0 as std::ffi::c_int as std::ffi::c_long
            && *T.offset((n - 1 as std::ffi::c_int as std::ffi::c_long) as isize)
                as std::ffi::c_int != 0 as std::ffi::c_int
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *PLCP
                .offset(0 as std::ffi::c_int as isize) = 0 as std::ffi::c_int as int64_t;
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    libsais16x64_compute_phi_omp(SA, PLCP, n, 1 as std::ffi::c_int as sa_sint_t);
    libsais16x64_compute_plcp_gsa_omp(T, PLCP, n, 1 as std::ffi::c_int as sa_sint_t);
    0 as std::ffi::c_int as int64_t
}
#[no_mangle]
pub unsafe extern "C" fn libsais16x64_lcp(
    mut PLCP: *const int64_t,
    mut SA: *const int64_t,
    mut LCP: *mut int64_t,
    mut n: int64_t,
) -> int64_t {
    if PLCP.is_null() || SA.is_null() || LCP.is_null()
        || n < 0 as std::ffi::c_int as std::ffi::c_long
    {
        return -(1 as std::ffi::c_int) as int64_t
    } else if n <= 1 as std::ffi::c_int as std::ffi::c_long {
        if n == 1 as std::ffi::c_int as std::ffi::c_long {
            *LCP
                .offset(
                    0 as std::ffi::c_int as isize,
                ) = *PLCP.offset(*SA.offset(0 as std::ffi::c_int as isize) as isize);
        }
        return 0 as std::ffi::c_int as int64_t;
    }
    libsais16x64_compute_lcp_omp(PLCP, SA, LCP, n, 1 as std::ffi::c_int as sa_sint_t);
    0 as std::ffi::c_int as int64_t
}
