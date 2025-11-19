#pragma once

// --- SIMPLE ---
// A basis defines the building blocks of the algebra
// Typically represented e1, e2, e3... eN
// These basis are our axis. Ex: x = e1, y = e2, z = e3
//
// We can combine basis to construct blades
// Blades are constructed by taking the outer product of one or more basis
// e1 (1-vector) - a vector along axis e1
// e2^e3 (2-vector) - an oriented plane along axis e2 e3
// e1^e2^e3 (3-vector) - an oriented cube along axis e1 e2 e3
// etc...
//
// Blades are represented using bitmasks