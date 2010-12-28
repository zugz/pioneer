#include "StarSystem.h"
#include "Sector.h"
#include "custom_starsystems.h"
#include "Serializer.h"
#include "NameGenerator.h"

#define CELSIUS	273.15
//#define DEBUG_DUMP

// minimum moon mass a little under Europa's
static const sfloat MIN_MOON_MASS = sfloat(1,30000); // earth masses
static const sfloat MIN_MOON_DIST = sfloat(15,10000); // AUs
static const sfloat MAX_MOON_DIST = sfloat(2, 100); // AUs
// if binary stars have separation s, planets can have stable
// orbits at (0.5 * s * SAFE_DIST_FROM_BINARY)
static const sfloat SAFE_DIST_FROM_BINARY = sfloat(5,1);
static const sfloat PLANET_MIN_SEPARATION = sfloat(135,100);

// very crudely
static const sfloat AU_SOL_RADIUS = sfloat(305,65536);
static const sfloat AU_EARTH_RADIUS = sfloat(3, 65536);

static const struct {
	sfloat meltingPoint;
	sfloat boilingPoint; // at 1 bar
	sfloat enthalpyVap;
	char irAbsorption[20]; // percent. 2 micrometer wide bands centred on 1um, 3um, 5um, ..., 39um
} s_chemStats[CHEM_MAX] = {
	/* H2 */
	{ sfloat(14,1), sfloat(20,1), sfloat(449,1) },
	/* O2 */
	{ sfloat(54,1), sfloat(90,1), sfloat(6820,1), // IR absorption actually for ozone
          { 0, 10, 7, 0, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
	/* N2 */
	{ sfloat(63,1), sfloat(77,1), sfloat(2792,1) },
	/* H2O */
	{ sfloat(27315,100), sfloat(37315,100), sfloat(40657,1),
          { 5, 10, 30, 80, 5, 0, 0, 0, 20, 25, 40, 50, 60, 85, 99, 99, 99, 99, 99, 99 } },
	/* CO2 */
	{ sfloat(195,1), sfloat(216,1), sfloat(15326,1),
          { 0, 15, 15, 0, 0, 0, 15, 99, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
	/* CH4 */
	{ sfloat(91,1), sfloat(112,1), sfloat(8180,1),
          { 0, 10, 0, 15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } },
	/* NH3 */
	{ sfloat(195,1), sfloat(240,1), sfloat(23350,1),

	}
};

// indexed by enum type turd  
float StarSystem::starColors[][3] = {
	{ 0, 0, 0 }, // gravpoint
	{ 0.5, 0.0, 0.0 }, // brown dwarf
	{ 1.0, 0.2, 0.0 }, // M
	{ 1.0, 0.6, 0.1 }, // K
	{ 0.4, 0.4, 0.8 }, // white dwarf
	{ 1.0, 1.0, 0.4 }, // G
	{ 1.0, 1.0, 0.8 }, // F
	{ 1.0, 1.0, 1.0 }, // A
	{ 0.7, 0.7, 1.0 }, // B
	{ 1.0, 0.7, 1.0 }  // O
};

// indexed by enum type turd  
float StarSystem::starRealColors[][3] = {
	{ 0, 0, 0 }, // gravpoint
	{ 0.5, 0.0, 0.0 }, // brown dwarf
	{ 1.0, 0.5, 0.2 }, // M
	{ 1.0, 1.0, 0.4 }, // K
	{ 1.0, 1.0, 1.0 }, // white dwarf
	{ 1.0, 1.0, 0.95 }, // G
	{ 1.0, 1.0, 1.0 }, // F
	{ 1.0, 1.0, 1.0 }, // A
	{ 0.8, 0.8, 1.0 }, // B
	{ 1.0, 0.8, 1.0 }  // O
};

float StarSystem::starLuminosities[] = {
	0,
	0.0003f, // brown dwarf
	0.08f, // M0
	0.38f, // K0
	1.2f, // G0
	5.1f, // F0
	24.0f, // A0
	24000.0f, // B0
	200000.0f, // O5
};

static const struct SBodySubTypeInfo {
	SBody::BodySuperType supertype;
	int mass[2]; // min,max % sol for stars, unused for planets
	int radius; // % sol radii for stars, % earth radii for planets
	const char *description;
	const char *icon;
	int tempMin, tempMax;
} bodyTypeInfo[SBody::TYPE_MAX] = {
	{
		SBody::SUPERTYPE_NONE, {}, 0, "Shouldn't see this!",
	}, {
		SBody::SUPERTYPE_STAR,
		{2,8}, 30, "Brown dwarf sub-stellar object",
		"icons/object_brown_dwarf.png",
		1000, 2000
	}, {
		SBody::SUPERTYPE_STAR,
		{10,47}, 50, "Type 'M' red star",
		"icons/object_star_m.png",
		2000, 3500
	}, {
		SBody::SUPERTYPE_STAR,
		{50,78}, 90, "Type 'K' orange star",
		"icons/object_star_k.png",
		3500, 5000
	}, {
		SBody::SUPERTYPE_STAR,
		{20,100}, 1, "White dwarf",
		"icons/object_white_dwarf.png",
		4000, 40000
	}, { 
		SBody::SUPERTYPE_STAR,
		{80,110}, 110, "Type 'G' yellow star",
		"icons/object_star_g.png",
		5000, 6000
	}, {
		SBody::SUPERTYPE_STAR,
		{115,170}, 140, "Type 'F' white star",
		"icons/object_star_f.png",
		6000, 7500
	}, {
		SBody::SUPERTYPE_STAR,
		{180,320}, 210, "Type 'A' hot white star",
		"icons/object_star_a.png",
		7500, 10000
	}, {
		SBody::SUPERTYPE_STAR,
		{400,1800}, 700, "Bright type 'B' blue star",
		"icons/object_star_b.png",
		10000, 30000
	}, {
		SBody::SUPERTYPE_STAR,
		{2000,4000}, 1600, "Hot, massive type 'O' blue star",
		"icons/object_star_o.png",
		30000, 60000
	}, {
		SBody::SUPERTYPE_GAS_GIANT,
		{}, 390, "Small gas giant",
		"icons/object_planet_small_gas_giant.png"
	}, {
		SBody::SUPERTYPE_GAS_GIANT,
		{}, 950, "Medium gas giant",
		"icons/object_planet_medium_gas_giant.png"
	}, {
		SBody::SUPERTYPE_GAS_GIANT,
		{}, 1110, "Large gas giant",
		"icons/object_planet_large_gas_giant.png"
	}, {
		SBody::SUPERTYPE_GAS_GIANT,
		{}, 1500, "Very large gas giant",
		"icons/object_planet_large_gas_giant.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 1, "Asteroid",
		"icons/object_planet_asteroid.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 2, "Large asteroid",
		"icons/object_planet_asteroid.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 26, "Small, rocky dwarf planet", // moon radius
		"icons/object_planet_dwarf.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 52, "Small, rocky planet with a thin atmosphere", // mars radius
		"icons/object_planet_small.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky frozen planet with a thin nitrogen atmosphere", // earth radius
		"icons/object_planet_water_n2.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Dead world that once housed it's own intricate ecosystem.", // earth radius
		"icons/object_planet_desert.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a carbon dioxide atmosphere",
		"icons/object_planet_co2.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a methane atmosphere",
		"icons/object_planet_methane.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Water world with vast oceans and a thick nitrogen atmosphere",
		"icons/object_planet_water_n1.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a thick carbon dioxide atmosphere",
		"icons/object_planet_co2.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Rocky planet with a thick methane atmosphere",
		"icons/object_planet_methane.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "Highly volcanic world",
		"icons/object_planet_volcanic.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 100, "World with indigenous life and an oxygen atmosphere",
		"icons/object_planet_life.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 60, "Marginal terraformed world with minimal plant life",
		"icons/object_planet_life3.png"
	}, {
		SBody::SUPERTYPE_ROCKY_PLANET,
		{}, 90, "Fully terraformed world with introduced species from numerous successful colonies",
		"icons/object_planet_life2.png"
	}, {
		SBody::SUPERTYPE_STARPORT,
		{}, 0, "Orbital starport",
		"icons/object_orbital_starport.png"
	}, {
		SBody::SUPERTYPE_STARPORT,
		{}, 0, "Starport",
	}
};

SBody::BodySuperType SBody::GetSuperType() const
{
	return bodyTypeInfo[type].supertype;
}

const char *SBody::GetAstroDescription()
{
	return bodyTypeInfo[type].description;
}

const char *SBody::GetIcon()
{
	return bodyTypeInfo[type].icon;
}

/*
 * Position a surface starport anywhere. Space.cpp::MakeFrameFor() ensures it
 * is on dry land (discarding this position if necessary)
 */
static void position_settlement_on_planet(SBody *b)
{
	MTRand r(b->seed);
	// used for orientation on planet surface
	b->orbit.rotMatrix = matrix4x4d::RotateZMatrix(2*M_PI*r.Double()) *
			matrix4x4d::RotateYMatrix(2*M_PI*r.Double());
}

double SBody::GetMaxChildOrbitalDistance() const
{
	double max = 0;
	for (unsigned int i=0; i<children.size(); i++) {
		if (children[i]->orbMax.ToDouble() > max) {
			max = children[i]->orbMax.ToDouble();	
		}
	}
	return AU * max;
}


static inline Sint64 isqrt(Sint64 a)
{
	Sint64 ret=0;
	Sint64 s;
	Sint64 ret_sq=-a-1;
	for(s=62; s>=0; s-=2){
		Sint64 b;
		ret+= ret;
		b=ret_sq + ((2*ret+1)<<s);
		if(b<0){
			ret_sq=b;
			ret++;
		}
	}
	return ret;
}


/*
 * These are the nice floating point surface temp calculating turds.
 *
static const double boltzman_const = 5.6704e-8;
static double calcEnergyPerUnitAreaAtDist(double star_radius, double star_temp, double object_dist)
{
	const double total_solar_emission = boltzman_const *
		star_temp*star_temp*star_temp*star_temp*
		4*M_PI*star_radius*star_radius;

	return total_solar_emission / (4*M_PI*object_dist*object_dist);
}

// bond albedo, not geometric
static double CalcSurfaceTemp(double star_radius, double star_temp, double object_dist, double albedo, double greenhouse)
{
	const double energy_per_meter2 = calcEnergyPerUnitAreaAtDist(star_radius, star_temp, object_dist);
	const double surface_temp = pow(energy_per_meter2*(1-albedo)/(4*(1-greenhouse)*boltzman_const), 0.25);
	return surface_temp;
}
*/
/*
 * Instead we use these butt-ugly overflow-prone spat of ejaculate:
 */
/*
 * star_radius in sol radii
 * star_temp in kelvin,
 * object_dist in AU
 * return Watts/m^2
 */
static sfloat calcEnergyPerUnitAreaAtDist(sfloat star_radius, int star_temp, sfloat object_dist)
{
	sfloat temp = star_temp * sfloat(1,10000);
	const sfloat total_solar_emission =
		temp*temp*temp*temp*star_radius*star_radius;
	
	return sfloat(1744665451,100000)*(total_solar_emission / (object_dist*object_dist));
}

static int CalcSurfaceTemp(const SBody *primary, sfloat distToPrimary, sfloat albedo, sfloat greenhouse)
{
	sfloat energy_per_meter2;
	if (primary->type == SBody::TYPE_GRAVPOINT) {
		// binary. take energies of both stars
		energy_per_meter2 = calcEnergyPerUnitAreaAtDist(primary->children[0]->radius,
			primary->children[0]->averageTemp, distToPrimary);
		energy_per_meter2 += calcEnergyPerUnitAreaAtDist(primary->children[1]->radius,
			primary->children[1]->averageTemp, distToPrimary);
	} else {
		energy_per_meter2 = calcEnergyPerUnitAreaAtDist(primary->radius, primary->averageTemp, distToPrimary);
	}
	const sfloat surface_temp_pow4 = energy_per_meter2*(1-albedo)/(1-greenhouse);
	return sfloat::Pow(surface_temp_pow4 * 4409673, sfloat(1,-2,false)).ToInt32();
}

vector3d Orbit::OrbitalPosAtTime(double t)
{
	const double e = eccentricity;
	// mean anomaly
	const double M = 2*M_PI*t / period;
	// eccentric anomaly
	// NR method to solve for E: M = E-sin(E)
	double E = M;
	for (int iter=5; iter > 0; --iter) {
		E = E - (E-e*(sin(E))-M) / (1.0 - e*cos(E));
	}
	// heliocentric distance
	double r = semiMajorAxis * (1.0 - e*cos(E));
	// true anomaly (angle of orbit position)
	double cos_v = (cos(E) - e) / (1.0 - e*cos(E));
	double sin_v = (sqrt(1.0-e*e)*sin(E))/ (1.0 - e*cos(E));

	vector3d pos = vector3d(-cos_v*r, sin_v*r, 0);
	pos = rotMatrix * pos;
	return pos;
}

vector3d Orbit::EvenSpacedPosAtTime(double t)
{
	const double e = eccentricity;
	const double M = 2*M_PI*t;
	const double v = 2*atan(sqrt((1+e)/(1-e)) * tan(M/2.0));
	const double r = semiMajorAxis * (1 - e*e) / (1 + e*cos(v));
	vector3d pos = vector3d(-cos(v)*r, sin(v)*r, 0);
	pos = rotMatrix * pos;
	return pos;
}

double calc_orbital_period(double semiMajorAxis, double centralMass)
{
	return 2.0*M_PI*sqrt((semiMajorAxis*semiMajorAxis*semiMajorAxis)/(G*centralMass));
}

EXPORT_OOLUA_FUNCTIONS_0_NON_CONST(SBodyPath)
EXPORT_OOLUA_FUNCTIONS_4_CONST(SBodyPath,
		GetBodyName, GetSeed, GetSystem, GetParent)

SBodyPath::SBodyPath(): SysLoc()
{
	sbodyId = 0;
}
SBodyPath::SBodyPath(int sectorX, int sectorY, int systemNum): SysLoc(sectorX, sectorY, systemNum)
{
	sbodyId = 0;
}

void SBodyPath::Serialize(Serializer::Writer &wr) const
{
	SysLoc::Serialize(wr);
	wr.Int32(sbodyId);
}

void SBodyPath::Unserialize(Serializer::Reader &rd, SBodyPath *path)
{
	SysLoc::Unserialize(rd, path);
	path->sbodyId = rd.Int32();
}

const char *SBodyPath::GetBodyName() const
{
	return GetSBody()->name.c_str();
}

Uint32 SBodyPath::GetSeed() const
{
	return GetSBody()->seed;
}

const SBody *SBodyPath::GetSBody() const
{
	StarSystem *s = StarSystem::GetCached(*this);
	return s->GetBodyByPath(this);
}
	
SBodyPath *SBodyPath::GetParent() const {
	SBodyPath *p = new SBodyPath;
	StarSystem *sys = StarSystem::GetCached(*this);
	sys->GetPathOf(sys->GetBodyByPath(this)->parent, p);
	return p;
}

template <class T>
static void shuffle_array(MTRand &rand, T *array, int len)
{
	for (int i=0; i<len; i++) {
		int pos = rand.Int32(len);
		T temp = array[i];
		array[i] = array[pos];
		array[pos] = temp;
	}
}

/*
 * Doesn't try very hard
 */
bool StarSystem::GetRandomStarportNearButNotIn(MTRand &rand, SBodyPath *outDest) const
{
	int sx = this->SectorX() + rand.Int32(3) - 1;
	int sy = this->SectorY() + rand.Int32(3) - 1;
	Sector sec(sx, sy);
	const int numSys = sec.m_systems.size();
	int *idxs = new int[numSys];
	// examine the systems in random order
	for (int i=0; i<numSys; i++) idxs[i] = i;
	shuffle_array<int>(rand, idxs, numSys);

	for (int i=0; i<numSys; i++) {
		if ((sx == this->SectorX()) &&
		    (sy == this->SectorY()) &&
		    (idxs[i] == this->SystemIdx())) continue;

		StarSystem *sys = new StarSystem(sx, sy, idxs[i]);

		const int numStations = sys->m_spaceStations.size();
		if (numStations) {
			sys->GetPathOf(sys->m_spaceStations[rand.Int32(numStations)],
					outDest);
			delete sys;
			return true;
		}
		delete sys;
	}
	return false;
}

SBody *StarSystem::GetBodyByPath(const SBodyPath *path) const
{
	assert((m_loc.sectorX == path->sectorX) || (m_loc.sectorY == path->sectorY) ||
	       (m_loc.systemNum == path->systemNum));
	assert(path->sbodyId < m_bodies.size());

	return m_bodies[path->sbodyId];
}

void StarSystem::GetPathOf(const SBody *sbody, SBodyPath *path) const
{
	*path = SBodyPath();

	path->sectorX = m_loc.sectorX;
	path->sectorY = m_loc.sectorY;
	path->systemNum = m_loc.systemNum;
	path->sbodyId = sbody->id;
}

/*
struct CustomSBody {
	const char *name; // null to end system
	SBody::BodyType type;
	int primaryIdx;  // -1 for primary
	sfloat radius; // in earth radii for planets, sol radii for stars
	sfloat mass; // earth masses or sol masses
	int averageTemp; // kelvin
	sfloat semiMajorAxis; // in AUs
	sfloat eccentricity;
};
*/
void StarSystem::CustomGetKidsOf(SBody *parent, const CustomSBody *customDef, const int primaryIdx, int *outHumanInfestedness, MTRand &rand)
{
	const CustomSBody *c = customDef;
	for (int i=0; c->name; c++, i++) {
		if (c->primaryIdx != primaryIdx) continue;
		
		SBody *kid = NewBody();
		SBody::BodyType type = c->type;
		kid->seed = rand.Int32();
		kid->type = type;
		kid->parent = parent;
		kid->radius = c->radius;
		kid->mass = c->mass;
		kid->averageTemp = c->averageTemp;
		kid->name = c->name;
		kid->rotationPeriod = c->rotationPeriod;
		kid->eccentricity = c->eccentricity;
		kid->axialTilt = c->axialTilt;
		kid->semiMajorAxis = c->semiMajorAxis;
		kid->orbit.eccentricity = c->eccentricity.ToDouble();
		kid->orbit.semiMajorAxis = c->semiMajorAxis.ToDouble() * AU;
		kid->orbit.period = calc_orbital_period(kid->orbit.semiMajorAxis, parent->GetMass());
		kid->heightMapFilename = c->heightMapFilename;

		if (kid->type == SBody::TYPE_STARPORT_SURFACE) {
			kid->orbit.rotMatrix = matrix4x4d::RotateYMatrix(c->longitude) *
				matrix4x4d::RotateXMatrix(-0.5*M_PI + c->latitude);
		} else {
			if (kid->orbit.semiMajorAxis < 1.2 * parent->GetRadius()) {
				Error("%s's orbit is too close to its parent", c->name);
			}
			kid->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rand.Double(2*M_PI)) *
				matrix4x4d::RotateXMatrix(-0.5*M_PI + c->latitude);
		}
		if (kid->GetSuperType() == SBody::SUPERTYPE_STARPORT) {
			(*outHumanInfestedness)++;
		}
		parent->children.push_back(kid);

		// perihelion and aphelion (in AUs)
		kid->orbMin = c->semiMajorAxis - c->eccentricity*c->semiMajorAxis;
		kid->orbMax = 2*c->semiMajorAxis - kid->orbMin;

		CustomGetKidsOf(kid, customDef, i, outHumanInfestedness, rand);
	}
}

void StarSystem::GenerateFromCustom(const CustomSystem *customSys, MTRand &rand)
{
	// find primary
	const CustomSBody *csbody = customSys->sbodies;

	int idx = 0;
	while ((csbody->name) && (csbody->primaryIdx != -1)) { csbody++; idx++; }
	assert(csbody->primaryIdx == -1);

	rootBody = NewBody();
	SBody::BodyType type = csbody->type;
	rootBody->type = type;
	rootBody->parent = NULL;
	rootBody->seed = rand.Int32();
	rootBody->radius = csbody->radius;
	rootBody->mass = csbody->mass;
	rootBody->averageTemp = csbody->averageTemp;
	rootBody->name = csbody->name;
	
	int humanInfestedness = 0;
	CustomGetKidsOf(rootBody, customSys->sbodies, idx, &humanInfestedness, rand);
	Populate(false);

}

void StarSystem::MakeStarOfType(SBody *sbody, SBody::BodyType type, MTRand &rand)
{
	sbody->type = type;
	sbody->radius = sfloat(bodyTypeInfo[type].radius, 100);
	sbody->mass = sfloat(rand.Int32(bodyTypeInfo[type].mass[0],
				bodyTypeInfo[type].mass[1]), 100);
	sbody->averageTemp = rand.Int32(bodyTypeInfo[type].tempMin,
				bodyTypeInfo[type].tempMax);
}

void StarSystem::MakeRandomStar(SBody *sbody, MTRand &rand)
{
	SBody::BodyType type = (SBody::BodyType)rand.Int32((int)SBody::TYPE_STAR_MIN, (int)SBody::TYPE_STAR_MAX);
	MakeStarOfType(sbody, type, rand);
}

void StarSystem::MakeStarOfTypeLighterThan(SBody *sbody, SBody::BodyType type, sfloat maxMass, MTRand &rand)
{
	int tries = 16;
	do {
		MakeStarOfType(sbody, type, rand);
	} while ((sbody->mass > maxMass) && (--tries));
}

void StarSystem::MakeBinaryPair(SBody *a, SBody *b, sfloat minDist, MTRand &rand)
{
	sfloat m = a->mass + b->mass;
	sfloat a0 = b->mass / m;
	sfloat a1 = a->mass / m;
	a->eccentricity = rand.NSfloat(3);
	int mul = 1;

	do {
		switch (rand.Int32(3)) {
			case 2: a->semiMajorAxis = sfloat(rand.Int32(100,10000), 100); break;
			case 1: a->semiMajorAxis = sfloat(rand.Int32(10,1000), 100); break;
			default:
			case 0: a->semiMajorAxis = sfloat(rand.Int32(1,100), 100); break;
		}
		a->semiMajorAxis *= mul;
		mul *= 2;
	} while (a->semiMajorAxis < minDist);

	a->orbit.eccentricity = a->eccentricity.ToDouble();
	a->orbit.semiMajorAxis = AU * (a->semiMajorAxis * a0).ToDouble();
	a->orbit.period = 60*60*24*365* a->semiMajorAxis.ToDouble() * sqrt(a->semiMajorAxis.ToDouble() / m.ToDouble());
	
	const float rotX = -0.5*M_PI;//(float)(rand.Double()*M_PI/2.0);
	const float rotY = (float)rand.Double(M_PI);
	a->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rotY) * matrix4x4d::RotateXMatrix(rotX);
	b->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rotY-M_PI) * matrix4x4d::RotateXMatrix(rotX);

	b->orbit.eccentricity = a->eccentricity.ToDouble();
	b->orbit.semiMajorAxis = AU * (a->semiMajorAxis * a1).ToDouble();
	b->orbit.period = a->orbit.period;
	
	sfloat orbMin = a->semiMajorAxis - a->eccentricity*a->semiMajorAxis;
	sfloat orbMax = 2*a->semiMajorAxis - orbMin;
	a->orbMin = orbMin;
	b->orbMin = orbMin;
	a->orbMax = orbMax;
	b->orbMax = orbMax;
}

SBody::SBody()
{
	heightMapFilename = 0;
}

/*
 * As my excellent comrades have pointed out, choices that depend on floating
 * point crap will result in different universes on different platforms.
 *
 * We must be sneaky and avoid floating point in these places.
 */
StarSystem::StarSystem(int sector_x, int sector_y, int system_idx)
{
	unsigned long _init[5] = { system_idx, sector_x, sector_y, UNIVERSE_SEED, 0 };
	memset(m_tradeLevel, 0, sizeof(m_tradeLevel));
	m_loc.sectorX = sector_x;
	m_loc.sectorY = sector_y;
	m_loc.systemNum = system_idx;
	rootBody = 0;
	if (system_idx == -1) return;

	Sector s = Sector(sector_x, sector_y);
	if (system_idx >= s.m_systems.size()) return;
	m_seed = s.m_systems[system_idx].seed;
	m_name = s.m_systems[system_idx].name;
	_init[4] = m_seed;
	MTRand rand;
	rand.seed(_init, 5);

	if (s.m_systems[system_idx].customSys) {
		const CustomSystem *custom = s.m_systems[system_idx].customSys;
		if (custom->shortDesc) m_shortDesc = custom->shortDesc;
		if (custom->longDesc) m_longDesc = custom->longDesc;
		if (custom->sbodies) {
			GenerateFromCustom(s.m_systems[system_idx].customSys, rand);
			return;
		}
	}

	SBody *star[4];
	SBody *centGrav1, *centGrav2;

	const int numStars = s.m_systems[system_idx].numStars;
	assert((numStars >= 1) && (numStars <= 4));

	if (numStars == 1) {
		SBody::BodyType type = s.m_systems[system_idx].starType[0];
		star[0] = NewBody();
		star[0]->parent = NULL;
		star[0]->name = s.m_systems[system_idx].name;
		star[0]->orbMin = 0;
		star[0]->orbMax = 0;
		MakeStarOfType(star[0], type, rand);
		rootBody = star[0];
		m_numStars = 1;
	} else {
		centGrav1 = NewBody();
		centGrav1->type = SBody::TYPE_GRAVPOINT;
		centGrav1->parent = NULL;
		centGrav1->name = s.m_systems[system_idx].name+" A,B";
		rootBody = centGrav1;

		SBody::BodyType type = s.m_systems[system_idx].starType[0];
		star[0] = NewBody();
		star[0]->name = s.m_systems[system_idx].name+" A";
		star[0]->parent = centGrav1;
		MakeStarOfType(star[0], type, rand);
		
		star[1] = NewBody();
		star[1]->name = s.m_systems[system_idx].name+" B";
		star[1]->parent = centGrav1;
		MakeStarOfTypeLighterThan(star[1], s.m_systems[system_idx].starType[1],
				star[0]->mass, rand);

		centGrav1->mass = star[0]->mass + star[1]->mass;
		centGrav1->children.push_back(star[0]);
		centGrav1->children.push_back(star[1]);
try_that_again_guvnah:
		MakeBinaryPair(star[0], star[1], sfloat(0), rand);

		m_numStars = 2;

		if (numStars > 2) {
			if (star[0]->orbMax > sfloat(100,1)) {
				// reduce to < 100 AU...
				goto try_that_again_guvnah;
			}
			// 3rd and maybe 4th star
			if (numStars == 3) {
				star[2] = NewBody();
				star[2]->name = s.m_systems[system_idx].name+" C";
				star[2]->orbMin = 0;
				star[2]->orbMax = 0;
				MakeStarOfTypeLighterThan(star[2], s.m_systems[system_idx].starType[2],
					star[0]->mass, rand);
				centGrav2 = star[2];
				m_numStars = 3;
			} else {
				centGrav2 = NewBody();
				centGrav2->type = SBody::TYPE_GRAVPOINT;
				centGrav2->name = s.m_systems[system_idx].name+" C,D";
				centGrav2->orbMax = 0;

				star[2] = NewBody();
				star[2]->name = s.m_systems[system_idx].name+" C";
				star[2]->parent = centGrav2;
				MakeStarOfTypeLighterThan(star[2], s.m_systems[system_idx].starType[2],
					star[0]->mass, rand);
				
				star[3] = NewBody();
				star[3]->name = s.m_systems[system_idx].name+" D";
				star[3]->parent = centGrav2;
				MakeStarOfTypeLighterThan(star[3], s.m_systems[system_idx].starType[3],
					star[2]->mass, rand);

				MakeBinaryPair(star[2], star[3], sfloat(0), rand);
				centGrav2->mass = star[2]->mass + star[3]->mass;
				centGrav2->children.push_back(star[2]);
				centGrav2->children.push_back(star[3]);
				m_numStars = 4;
			}
			SBody *superCentGrav = NewBody();
			superCentGrav->type = SBody::TYPE_GRAVPOINT;
			superCentGrav->parent = NULL;
			superCentGrav->name = s.m_systems[system_idx].name;
			centGrav1->parent = superCentGrav;
			centGrav2->parent = superCentGrav;
			rootBody = superCentGrav;
			const sfloat minDist = star[0]->orbMax + star[2]->orbMax;
			MakeBinaryPair(centGrav1, centGrav2, 4*minDist, rand);
			superCentGrav->children.push_back(centGrav1);
			superCentGrav->children.push_back(centGrav2);

		}
	}

	for (int i=0; i<m_numStars; i++) MakePlanetsAround(star[i], rand);

	if (m_numStars > 1) MakePlanetsAround(centGrav1, rand);
	if (m_numStars == 4) MakePlanetsAround(centGrav2, rand);

	Populate(true);

#ifdef DEBUG_DUMP
	Dump();
#endif /* DEBUG_DUMP */
}

#ifdef DEBUG_DUMP
struct thing_t {
	SBody* obj;
	vector3d pos;
	vector3d vel;
};
void StarSystem::Dump()
{
	std::vector<SBody*> obj_stack;
	std::vector<vector3d> pos_stack;
	std::vector<thing_t> output;
	
	SBody *obj = rootBody;
	vector3d pos = vector3d(0.0);

	while (obj) {
		vector3d p2 = pos;
		if (obj->parent) {
			p2 = pos + obj->orbit.OrbitalPosAtTime(1.0);
			pos = pos + obj->orbit.OrbitalPosAtTime(0.0);
		}

		if ((obj->type != SBody::TYPE_GRAVPOINT) &&
		    (obj->GetSuperType() != SBody::SUPERTYPE_STARPORT)) {
			struct thing_t t;
			t.obj = obj;
			t.pos = pos;
			t.vel = (p2-pos);
			output.push_back(t);
		}
		for (std::vector<SBody*>::iterator i = obj->children.begin();
				i != obj->children.end(); ++i) {
			obj_stack.push_back(*i);
			pos_stack.push_back(pos);
		}
		if (obj_stack.size() == 0) break;
		pos = pos_stack.back();
		obj = obj_stack.back();
		pos_stack.pop_back();
		obj_stack.pop_back();
	}

	FILE *f = fopen("starsystem.dump", "w");
	fprintf(f, "%d bodies\n", output.size());
	fprintf(f, "0 steps\n");
	for (std::vector<thing_t>::iterator i = output.begin();
			i != output.end(); ++i) {
		fprintf(f, "B:%lf,%lf:%lf,%lf,%lf,%lf:%lf:%d:%lf,%lf,%lf\n",
				(*i).pos.x, (*i).pos.y, (*i).pos.z,
				(*i).vel.x, (*i).vel.y, (*i).vel.z,
				(*i).obj->GetMass(), 0,
				1.0, 1.0, 1.0);
	}
	fclose(f);
	printf("Junk dumped to starsystem.dump\n");
}
#endif /* DEBUG_DUMP */

/*
 * http://en.wikipedia.org/wiki/Hill_sphere
 */
sfloat SBody::CalcHillRadius() const
{
	if (GetSuperType() <= SUPERTYPE_STAR) {
		return sfloat(0);
	} else {
		// masses in earth masses
		return semiMajorAxis * (1 - eccentricity) * sfloat::CubeRoot(mass / (3 * parent->GetMassInEarths()));
	}
}

static sfloat mass_from_disk_area(sfloat a, sfloat b, sfloat max)
{
	// so, density of the disk with distance from star goes like so: 1 - x/discMax
	//
	// --- 
	//    ---
	//       --- <- zero at discMax
	//
	// Which turned into a disc becomes 2*pi*x - (2*pi*x*x)/discMax
	// Integral of which is: pi*x*x - (2/(3*discMax))*pi*x*x*x
	//
	// Because get_disc_density divides total_mass by
	// mass_from_disk_area(0, discMax, discMax) to find density, the
	// constant factors (pi) in this equation drop out.
	//
	b = (b > max ? max : b);
	assert(b>=a);
	assert(a<=max);
	assert(b<=max);
	assert(a>=0);
	sfloat one_over_3max = sfloat(2,1)/(3*max);
	return (b*b - one_over_3max*b*b*b) -
		(a*a - one_over_3max*a*a*a);
}

static sfloat get_disc_density(SBody *primary, sfloat discMin, sfloat discMax, sfloat percentOfPrimaryMass)
{
	discMax = MAX(discMax, discMin);
	sfloat total = mass_from_disk_area(discMin, discMax, discMax);
	return primary->GetMassInEarths() * percentOfPrimaryMass / total;
}

void StarSystem::MakePlanetsAround(SBody *primary, MTRand &rand)
{
	sfloat discMin = sfloat(0);
	sfloat discMax = sfloat(5000,1);
	sfloat discDensity;

	SBody::BodySuperType superType = primary->GetSuperType();

	if (superType <= SBody::SUPERTYPE_STAR) {
		if (primary->type == SBody::TYPE_GRAVPOINT) {
			/* around a binary */
			discMin = primary->children[0]->orbMax * SAFE_DIST_FROM_BINARY;
		} else {
			/* correct thing is roche limit, but lets ignore that because
			 * it depends on body densities and gives some strange results */
			discMin = 4 * primary->radius * AU_SOL_RADIUS;
		}
		if (primary->type == SBody::TYPE_WHITE_DWARF) {
			// white dwarfs will have started as stars < 8 solar
			// masses or so, so pick discMax according to that
			discMax = 100 * rand.NSfloat(2)*sfloat::Sqrt(sfloat(1,2) + sfloat(8,1)*rand.Sfloat());
		} else {
			discMax = 100 * rand.NSfloat(2)*sfloat::Sqrt(primary->mass);
		}
		// having limited discMin by bin-separation/fake roche, and
		// discMax by some relation to star mass, we can now compute
		// disc density
		discDensity = rand.Sfloat() * get_disc_density(primary, discMin, discMax, sfloat(2,100));

		if ((superType == SBody::SUPERTYPE_STAR) && (primary->parent)) {
			// limit planets out to 10% distance to star's binary companion
			discMax = MIN(discMax, primary->orbMin * sfloat(1,10));
		}

		/* in trinary and quaternary systems don't bump into other pair... */
		if (m_numStars >= 3) {
			discMax = MIN(discMax, sfloat(5,100)*rootBody->children[0]->orbMin);
		}
	} else {
		sfloat primary_rad = primary->radius * AU_EARTH_RADIUS;
		discMin = 4 * primary_rad;
		/* use hill radius to find max size of moon system. for stars botch it */
		discMax = MIN(discMax, sfloat(1,20)*primary->CalcHillRadius());
		
		discDensity = rand.Sfloat() * get_disc_density(primary, discMin, discMax, sfloat(1,500));
	}

	//sfloat discDensity = 20*rand.NSfloat(4);

	//printf("Around %s: Range %f -> %f AU\n", primary->name.c_str(), discMin.ToDouble(), discMax.ToDouble());

	sfloat initialJump = rand.NSfloat(5);
	sfloat pos = (sfloat(1,1) - initialJump)*discMin + (initialJump*discMax);

	while (pos < discMax) {
		// periapsis, apoapsis = closest, farthest distance in orbit
		sfloat periapsis = pos + pos*0.5*rand.NSfloat(2);/* + jump */;
		sfloat ecc = rand.NSfloat(3);
		sfloat semiMajorAxis = periapsis / (sfloat(1,1) - ecc);
		sfloat apoapsis = 2*semiMajorAxis - periapsis;
		if (apoapsis > discMax) break;

		sfloat mass;
		{
			const sfloat a = pos;
			const sfloat b = sfloat(135,100)*apoapsis;
			mass = mass_from_disk_area(a, b, discMax);
			mass *= rand.Sfloat() * discDensity;
		}

		SBody *planet = NewBody();
		planet->eccentricity = ecc;
		planet->axialTilt = sfloat(100,157)*rand.NSfloat(2);
		planet->semiMajorAxis = semiMajorAxis;
		planet->type = SBody::TYPE_PLANET_DWARF;
		planet->seed = rand.Int32();
		planet->tmp = 0;
		planet->parent = primary;
	//	planet->radius = EARTH_RADIUS*bodyTypeInfo[type].radius;
		planet->mass = mass;
		planet->rotationPeriod = sfloat(rand.Int32(1,200), 24);

		planet->orbit.eccentricity = ecc.ToDouble();
		planet->orbit.semiMajorAxis = semiMajorAxis.ToDouble() * AU;
		planet->orbit.period = calc_orbital_period(planet->orbit.semiMajorAxis, primary->GetMass());
		planet->orbit.rotMatrix = matrix4x4d::RotateYMatrix(rand.Double(2*M_PI)) *
			matrix4x4d::RotateXMatrix(-0.5*M_PI + rand.NDouble(5)*M_PI/2.0);
		planet->orbMin = periapsis;
		planet->orbMax = apoapsis;
		primary->children.push_back(planet);

		/* minimum separation between planets of 1.35 */
		pos = apoapsis * sfloat(135,100);
	}

	int idx=0;
	bool make_moons = superType <= SBody::SUPERTYPE_STAR;
	
	for (std::vector<SBody*>::iterator i = primary->children.begin(); i != primary->children.end(); ++i) {
		// planets around a binary pair [gravpoint] -- ignore the stars...
		if ((*i)->GetSuperType() == SBody::SUPERTYPE_STAR) continue;
		// Turn them into something!!!!!!!
		char buf[8];
		if (superType <= SBody::SUPERTYPE_STAR) {
			// planet naming scheme
			snprintf(buf, sizeof(buf), " %c", 'b'+idx);
		} else {
			// moon naming scheme
			snprintf(buf, sizeof(buf), " %d", 1+idx);
		}
		(*i)->name = primary->name+buf;
		(*i)->PickPlanetType(this, rand);
		if (make_moons) MakePlanetsAround(*i, rand);
		idx++;
	}
}

/*
 * For moons distance from star is not orbMin, orbMax.
 */
const SBody *SBody::FindStarAndTrueOrbitalRange(sfloat &orbMin, sfloat &orbMax)
{
	const SBody *planet = this;
	const SBody *star = this->parent;

	assert(star);

	/* while not found star yet.. */
	while (star->GetSuperType() > SBody::SUPERTYPE_STAR) {
		planet = star;
		star = star->parent;
	}

	orbMin = planet->orbMin;
	orbMax = planet->orbMax;
	return star;
}

sfloat SBody::CalcGlobalWarming() const
{
	// emission by planet
	sfloat irEmission[20];
	// planck blackbody radiation formula:
	// I(l) = (2*h*c**2) / (l**5 * (e**((hc)/(lkT)) - 1))
	// Which without constants is:
	// I(l) = 1 / (l**5 * (e**(1/(lT))-1))
	// I(l) = power, l = wavelength, T = temp (kelvin)
	// l = 0.000069 corresponds to 1um
	const sfloat oneMicron = sfloat(69,1000000);
	const sfloat T = sfloat(averageTemp,1);

	printf("Temperature: %f K (%d)\n", T.ToDouble(), averageTemp);
	sfloat total;
	for (int i=0; i<20; i++) {
		const sfloat l = sfloat(2*i+1,1) * oneMicron;
		irEmission[i] = sfloat(1,1) / (l*l*l*l*l * (sfloat::Exp(sfloat(1,1) / (l*T)) - sfloat(1,1)));
		total += irEmission[i];
	}
	/* sortof normalize the shit */
	total = sfloat(1,1) / total;
	for (int i=0; i<20; i++) {
		printf("%.0f um emission: %f%%\n", (double)(2*i+1), 100.0*(total * irEmission[i]).ToDouble());
	}

	sfloat greenhouse;

	for (int band=0; band<20; band++) {
		sfloat bandabs;
		for (int gas=0; gas<CHEM_MAX; gas++) {
			bandabs += m_gases[gas] * sfloat(s_chemStats[gas].irAbsorption[band], 100);
		}
		bandabs = bandabs / (sfloat(1,1) + bandabs);
		printf("%.0f um absorption: %f%%\n", (double)(2*band+1), 100.0*(bandabs).ToDouble());
		greenhouse += total * irEmission[band] * bandabs;
	}
	printf("Greenhouse: %f\n", greenhouse.ToDouble());
	return sfloat((greenhouse*sfloat(1000,1)).ToInt32(), 1000);
}

void SBody::PickPlanetType(StarSystem *system, MTRand &rand)
{
#if 0
	///////////////////////
	// OK, pick composition
	///////////////////////
	crustCarbon = rand.NSfloat(2);
	crustOxygen = rand.Sfloat();
	crustSilicon = rand.Sfloat();
	crustLightEarths = rand.Sfloat();
	crustLightTrans = rand.NSfloat(2);
	crustHeavyMetals = rand.NSfloat(3);
	{
		sfloat invTotal = sfloat(1,1) / (crustCarbon + crustOxygen + crustSilicon +
			crustLightEarths + crustLightTrans + crustHeavyMetals);
		if (invTotal == 0) invTotal = sfloat(1,1);
		crustCarbon *= invTotal;
		crustOxygen *= invTotal;
		crustSilicon *= invTotal;
		crustLightEarths *= invTotal;
		crustLightTrans *= invTotal;
		crustHeavyMetals *= invTotal;
	}
	printf("============\n");
	printf("Crust: C %f, O %f, Si %f, Light Earths %f, Light Trans %f, Heavy metals %f\n",
			crustCarbon.ToDouble(),
			crustOxygen.ToDouble(),
			crustSilicon.ToDouble(),
			crustLightEarths.ToDouble(),
			crustLightTrans.ToDouble(),
			crustHeavyMetals.ToDouble());
	// make terrible estimate of radius
	radius = sfloat::CubeRoot(mass);

	for (int i=0; i<CHEM_MAX; i++) {
		m_gases[i] = m_liquids[i] = m_ices[i] = 0;
	}

	{
		// inputs
		sfloat carbon = rand.NSfloat(3);
		m_gases[CHEM_H2] = rand.Sfloat();
		m_gases[CHEM_O2] = rand.Sfloat();
		m_gases[CHEM_N2] = rand.Sfloat();
		// 'react' them...
#define REACT(outProduct, amount_a, molarMass_a, amount_b, molarMass_b) { \
		outProduct = MIN(sfloat(molarMass_a + molarMass_b, molarMass_a) * amount_a, \
				 sfloat(molarMass_a + molarMass_b, molarMass_b) * amount_b); \
		amount_a -= sfloat(molarMass_a, molarMass_a + molarMass_b) * outProduct; \
		amount_b -= sfloat(molarMass_b, molarMass_a + molarMass_b) * outProduct; }

		// my chemistry is a bit rusty, but perhaps this is the right order
		// for things to react....
		REACT(m_gases[CHEM_H2O], m_gases[CHEM_O2], 16, m_gases[CHEM_H2], 2);
		REACT(m_gases[CHEM_CO2], m_gases[CHEM_O2], 32, carbon, 12);
		REACT(m_gases[CHEM_CH4], carbon, 12, m_gases[CHEM_H2], 4);
		REACT(m_gases[CHEM_NH3], m_gases[CHEM_N2], 14, m_gases[CHEM_H2], 3);

		// this equation is kosher. It is newton's gravity over the surface area:
		// pressure = planet_mass * atmosphere_mass / (radius^4)
		// pressure in earth atmospheres,
		// planet_mass in earth masses,
		// atmosphere_mass in earth atmosphere masses
		// radius in earth radii
		sfloat minDistToStar, maxDistToStar, averageDistToStar;
		const SBody *star = FindStarAndTrueOrbitalRange(minDistToStar, maxDistToStar);
		averageDistToStar = (minDistToStar+maxDistToStar)*sfloat(1,2);
		printf("Dist to star %.3fAU, star temperature: %dK\n", averageDistToStar.ToDouble(), star->averageTemp);

		for (int iteration=0; iteration<20; iteration++) {
			// OK, so when you condense the atmosphere out the pressure falls and
			// boiling points change, so it must be done in iterations
			sfloat m = sfloat(1,10); 
			sfloat atmosphereMass = sfloat(0);
			for (int i=0; i<CHEM_MAX; i++) atmosphereMass += m_gases[i];
			m_surfacePressure = (mass * atmosphereMass) / (radius*radius*radius*radius);

			sfloat albedo = sfloat(14,100);
			// H2O changes albedo by cloud formation
			albedo += sfloat::Sqrt(m_gases[CHEM_H2O]);
			// all ices raise albedo
			for (int i=0; i<CHEM_MAX; i++) { albedo += m_ices[i]; }
			// squish range in a highly dubious made up way...
			albedo = albedo / (sfloat(1,1) + albedo);

			averageTemp = CalcSurfaceTemp(star, averageDistToStar, albedo, sfloat(0));
			// so calculating the global warming factor is sortof black magic and fraud.
			// It uses some crude weightings and a surface area distribution
			// and a rather suspect compressing function to ensure 0.0-1.0 range.
			sfloat globalwarming = this->CalcGlobalWarming();
				//(m_gases[CHEM_CO2] +
				  //            sfloat(228,10)*m_gases[CHEM_CH4] +
		//			      sfloat(7,1000)*m_gases[CHEM_H2O]) / (radius*radius);
		//	globalwarming = sfloat::SqrtOf(globalwarming) / (sfloat(2,100) + sfloat::SqrtOf(globalwarming));
			printf("Global warming: %f, albedo: %f, surface atmospheric pressure %f\n", globalwarming.ToDouble(),
					albedo.ToDouble(), m_surfacePressure.ToDouble());

			printf("Global warming: %f\n", globalwarming.ToDouble());
			averageTemp = CalcSurfaceTemp(star, averageDistToStar, albedo, globalwarming);
			printf("New average temp %d K\n", averageTemp);
			printf("Average surface temp %dK, pressure %f bar\n", averageTemp, m_surfacePressure.ToDouble());
			printf("State	H2	O2	N2	H2O	CO2	CH4	NH3");
			printf("\nGas");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_gases[i].ToDouble()); }
			printf("\nLiquid");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_liquids[i].ToDouble()); }
			printf("\nIces");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_ices[i].ToDouble()); }
			printf("\n");
			//printf("Stage %.3f, surface pressure %f atmospheres, temp %dK\n", m.ToDouble(), m_surfacePressure.ToDouble(), averageTemp);
			for (int i=0; i<CHEM_MAX; i++) {
				if (s_chemStats[i].meltingPoint > sfloat(averageTemp,1)) {
					// move to ices state
					m_ices[i] += m * (m_gases[i] + m_liquids[i]);
					m_gases[i] -= m * m_gases[i];
					m_liquids[i] -= m * m_liquids[i];
				} else {
					sfloat boilingPoint = sfloat(1,1) / (sfloat::Log(sfloat(1,1)/MAX(sfloat(1,100),m_surfacePressure)) * 
						(sfloat(8314,1000) / s_chemStats[i].enthalpyVap) +
						(sfloat(1,1)/s_chemStats[i].boilingPoint));
					//printf("%d boils at %f K on this joint\n", i, boilingPoint.ToDouble());
					if (boilingPoint > sfloat(averageTemp,1)) {
						// move to liquids state
						m_liquids[i] += m * (m_gases[i] + m_ices[i]);
						m_gases[i] -= m * m_gases[i];
						m_ices[i] -= m * m_ices[i];
					} else {
						// move to gases state
						m_gases[i] += m * (m_liquids[i] + m_ices[i]);
						m_liquids[i] -= m * m_liquids[i];
						m_ices[i] -= m * m_ices[i];
					}

				}
			}
		}
			printf("Average surface temp %dK, pressure %f bar\n", averageTemp, m_surfacePressure.ToDouble());
			printf("State	H2	O2	N2	H2O	CO2	CH4	NH3");
			printf("\nGas");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_gases[i].ToDouble()); }
			printf("\nLiquid");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_liquids[i].ToDouble()); }
			printf("\nIces");
			for (int i=0; i<CHEM_MAX; i++) { printf("\t%.3f", m_ices[i].ToDouble()); }
			printf("\n");
	}
#endif

	////////////////////////

	sfloat albedo = rand.Sfloat() * sfloat(1,2);
	sfloat globalwarming = rand.Sfloat() * sfloat(9,10);
	// light planets have bugger all atmosphere
	if (mass < 1) globalwarming *= mass;
	// big planets get high global warming due to thick atmos
	if (mass > 3) globalwarming *= (mass-2);
	globalwarming = CLAMP(globalwarming, sfloat(0), sfloat(95,100));

	sfloat minDistToStar, maxDistToStar, averageDistToStar;
	const SBody *star = FindStarAndTrueOrbitalRange(minDistToStar, maxDistToStar);
	averageDistToStar = (minDistToStar+maxDistToStar)*sfloat(1,2);

	/* this is all of course a total fucking joke and un-physical */
	int bbody_temp;
	bool fiddle = false;
	for (int i=0; i<10; i++) {
		bbody_temp = CalcSurfaceTemp(star, averageDistToStar, albedo, globalwarming);
		//printf("temp %f, albedo %f, globalwarming %f\n", bbody_temp, albedo, globalwarming);
		// extreme high temperature and low mass causes atmosphere loss
#define ATMOS_LOSS_MASS_CUTOFF	2
#define ATMOS_TEMP_CUTOFF	400
#define FREEZE_TEMP_CUTOFF	220
		if ((bbody_temp > ATMOS_TEMP_CUTOFF) &&
		   (mass < ATMOS_LOSS_MASS_CUTOFF)) {
		//	printf("atmos loss\n");
			globalwarming = globalwarming * (mass/ATMOS_LOSS_MASS_CUTOFF);
			fiddle = true;
		}
		if (!fiddle) break;
		fiddle = false;
	}
	// this is utter rubbish. should decide atmosphere composition and then freeze out
	// components of it in the previous loop
	if ((bbody_temp < FREEZE_TEMP_CUTOFF) && (mass < 5)) {
		globalwarming *= 0.2;
		albedo = rand.Sfloat()*sfloat(5,100) + 0.9;
	}
	bbody_temp = CalcSurfaceTemp(star, averageDistToStar, albedo, globalwarming);
//	printf("= temp %f, albedo %f, globalwarming %f\n", bbody_temp, albedo, globalwarming);

	averageTemp = bbody_temp;

	if (mass > 317*13) {
		// more than 13 jupiter masses can fuse deuterium - is a brown dwarf
		type = SBody::TYPE_BROWN_DWARF;
		// prevent mass exceeding 65 jupiter masses or so, when it becomes a star
		// XXX since TYPE_BROWN_DWARF is supertype star, mass is now in
		// solar masses. what a fucking mess
		mass = MIN(mass, sfloat(317*65, 1)) / 332998;
	} else if (mass > 300) {
		type = SBody::TYPE_PLANET_LARGE_GAS_GIANT;
	} else if (mass > 90) {
		type = SBody::TYPE_PLANET_MEDIUM_GAS_GIANT;
	} else if (mass > 6) {
		type = SBody::TYPE_PLANET_SMALL_GAS_GIANT;
	} else {
		// terrestrial planets
		if (mass < sfloat(1,20000)) {
			type = SBody::TYPE_PLANET_ASTEROID;
		} else if (mass < sfloat(1, 15000)) {
			type = SBody::TYPE_PLANET_LARGE_ASTEROID;
		} else if (mass < sfloat(2,1000)) {
			type = SBody::TYPE_PLANET_DWARF;
		} else if ((mass < sfloat(2,10)) && (globalwarming < sfloat(5,100))) {
			type = SBody::TYPE_PLANET_SMALL;
		} else if (mass < 3) {
			if ((averageTemp > CELSIUS-60) && (averageTemp < CELSIUS+200)) {
				// try for life
				int minTemp = CalcSurfaceTemp(star, maxDistToStar, albedo, globalwarming);
				int maxTemp = CalcSurfaceTemp(star, minDistToStar, albedo, globalwarming);

				if ((star->type != TYPE_BROWN_DWARF) &&
				    (star->type != TYPE_WHITE_DWARF) &&
				    (star->type != TYPE_STAR_O) &&
				    (minTemp > CELSIUS-10) && (minTemp < CELSIUS+60) &&
				    (maxTemp > CELSIUS-10) && (maxTemp < CELSIUS+60)) {
					type = SBody::TYPE_PLANET_INDIGENOUS_LIFE;
					humanActivity *= 4;
				} else if 
						 ((minTemp > CELSIUS-15) && (minTemp < CELSIUS+65) &&
						 (maxTemp > CELSIUS-15) && (maxTemp < CELSIUS+65)) {
						 type = SBody::TYPE_PLANET_TERRAFORMED_GOOD;
						 humanActivity *= 3;
						 }
				else if
						 ((minTemp > CELSIUS-25) && (minTemp < CELSIUS+70) &&
						 (maxTemp > CELSIUS-25) && (maxTemp < CELSIUS+70)) {
						 type = SBody::TYPE_PLANET_TERRAFORMED_POOR;
						 humanActivity *= 1;
						 }
				else if	
						((minTemp > CELSIUS-00) && (minTemp < CELSIUS+95) &&
						 (maxTemp > CELSIUS-00) && (maxTemp < CELSIUS+95)) {
						 type = SBody::TYPE_PLANET_WATER_THICK_ATMOS;
						 humanActivity *= 2;
						 }
				else if	
						((minTemp > CELSIUS-100) && (minTemp < CELSIUS+00) &&
						 (maxTemp > CELSIUS-100) && (maxTemp < CELSIUS+00)) {
						 type = SBody::TYPE_PLANET_WATER;
						 humanActivity *= 1;
						 }
				else if	
						((minTemp > CELSIUS+30) && (minTemp < CELSIUS+120) &&
						 (maxTemp > CELSIUS+30) && (maxTemp < CELSIUS+120)) {
						 type = SBody::TYPE_PLANET_DESERT;
						 humanActivity *= 1;
						//bit crap but it has the desired effect 
						}
				else
					{
					if (averageTemp > CELSIUS+10)
					{
					type = SBody::TYPE_PLANET_DESERT;
				    }
					else if (averageTemp < CELSIUS-10)
					{
					type = SBody::TYPE_PLANET_WATER;
					}}
			} else {
				if (rand.Int32(0,1)) type = SBody::TYPE_PLANET_CO2;
				else type = SBody::TYPE_PLANET_METHANE;
			}
		} else /* 3 < mass < 6 */ {
			if ((averageTemp > CELSIUS-5) && (averageTemp < CELSIUS+80)) {
				type = SBody::TYPE_PLANET_WATER_THICK_ATMOS;
			}
			else if ((averageTemp > CELSIUS-150) && (averageTemp < CELSIUS-3)) {
					type = SBody::TYPE_PLANET_WATER;
			}
			else if ((averageTemp > CELSIUS+70) && (averageTemp < CELSIUS+600)) {
					type = SBody::TYPE_PLANET_DESERT;
			}
			else {
				if (rand.Int32(0,1)) type = SBody::TYPE_PLANET_CO2_THICK_ATMOS;
				else type = SBody::TYPE_PLANET_METHANE_THICK_ATMOS;
			}
		}
		// kind of crappy
		if ((mass > sfloat(8,10)) && (!rand.Int32(0,15))) type = SBody::TYPE_PLANET_HIGHLY_VOLCANIC;
	}
	radius = sfloat(bodyTypeInfo[type].radius, 100);
}

void StarSystem::MakeShortDescription(MTRand &rand)
{
	m_econType = 0;
	if ((m_industrial > m_metallicity) && (m_industrial > m_agricultural)) {
		m_econType = ECON_INDUSTRY;
	} else if (m_metallicity > m_agricultural) {
		m_econType = ECON_MINING;
	} else {
		m_econType = ECON_AGRICULTURE;
	}

	/* Total population is in billions */
	if (m_totalPop == 0) {
		int dist = isqrt(1 + m_loc.sectorX*m_loc.sectorX + m_loc.sectorY*m_loc.sectorY);
		if (rand.Int32(dist) > 20) {
			m_shortDesc = "Unexplored system.";
		} else {
			m_shortDesc = "Small-scale prospecting. No registered settlements.";
		}
	} else if (m_totalPop < sfloat(1,10)) {
		switch (m_econType) {
			case ECON_INDUSTRY: m_shortDesc = "Small industrial outpost."; break;
			case ECON_MINING: m_shortDesc = "Some established mining."; break;
			case ECON_AGRICULTURE: m_shortDesc = "Young farming colony."; break;
		}
	} else if (m_totalPop < sfloat(1,2)) {
		switch (m_econType) {
			case ECON_INDUSTRY: m_shortDesc = "Industrial colony."; break;
			case ECON_MINING: m_shortDesc = "Mining colony."; break;
			case ECON_AGRICULTURE: m_shortDesc = "Outdoor agricultural world."; break;
		}
	} else if (m_totalPop < sfloat(5,1)) {
		switch (m_econType) {
			case ECON_INDUSTRY: m_shortDesc = "Heavy industry."; break;
			case ECON_MINING: m_shortDesc = "Extensive mining operations."; break;
			case ECON_AGRICULTURE: m_shortDesc = "Thriving outdoor world."; break;
		}
	} else {
		switch (m_econType) {
			case ECON_INDUSTRY: m_shortDesc = "Industrial hub system."; break;
			case ECON_MINING: m_shortDesc = "Vast strip-mining colony."; break;
			case ECON_AGRICULTURE: m_shortDesc = "High population outdoor world."; break;
		}
	}
}

/* percent */
#define MAX_COMMODITY_BASE_PRICE_ADJUSTMENT 25

void StarSystem::Populate(bool addSpaceStations)
{
	unsigned long _init[5] = { m_loc.systemNum, m_loc.sectorX, m_loc.sectorY, UNIVERSE_SEED };
	MTRand rand;
	rand.seed(_init, 4);

	/* Various system-wide characteristics */
	m_humanProx = sfloat(3,1) / isqrt(9 + 10*(m_loc.sectorX*m_loc.sectorX + m_loc.sectorY*m_loc.sectorY));
	m_metallicity = rand.Sfloat();
	m_techlevel = (m_humanProx*5).ToInt32() + rand.Int32(-2,2);
	m_techlevel = CLAMP(m_techlevel, 1, 5);
	m_econType = ECON_INDUSTRY;
	m_industrial = rand.Sfloat();
	m_agricultural = 0;

	/* system attributes */
	m_totalPop = sfloat(0);
	rootBody->PopulateStage1(this, m_totalPop);
	if (m_totalPop == 0) m_techlevel = 0;
	
//	printf("Trading rates:\n");
	// So now we have balances of trade of various commodities.
	// Lets use black magic to turn these into percentage base price
	// alterations
	int maximum = 0;
	for (int i=(int)Equip::FIRST_COMMODITY; i<=(int)Equip::LAST_COMMODITY; i++) {
		maximum = MAX(abs(m_tradeLevel[i]), maximum);
	}
	if (maximum) for (int i=(int)Equip::FIRST_COMMODITY; i<=(int)Equip::LAST_COMMODITY; i++) {
		m_tradeLevel[i] = (m_tradeLevel[i] * MAX_COMMODITY_BASE_PRICE_ADJUSTMENT) / maximum;
		m_tradeLevel[i] += rand.Int32(-5, 5);
	}
	
	for (int i=(int)Equip::FIRST_COMMODITY; i<=(int)Equip::LAST_COMMODITY; i++) {
		Equip::Type t = (Equip::Type)i;
		const EquipType &type = EquipType::types[t];
//		printf("%s: %d%%\n", type.name, m_tradeLevel[t]);
	}
//	printf("System total population %.3f billion, tech level %d\n", m_totalPop.ToFloat(), m_techlevel);
	Polit::GetSysPolitStarSystem(this, m_totalPop, m_polit);

	if (addSpaceStations) {
		rootBody->PopulateAddStations(this);
	}

	MakeShortDescription(rand);
}

/*
 * Set natural resources, tech level, industry strengths and population levels
 */
void SBody::PopulateStage1(StarSystem *system, sfloat &outTotalPop)
{
	for (unsigned int i=0; i<children.size(); i++) {
		children[i]->PopulateStage1(system, outTotalPop);
	}
	unsigned long _init[5] = { system->m_loc.systemNum, system->m_loc.sectorX,
			system->m_loc.sectorY, UNIVERSE_SEED, this->seed };
	MTRand rand;
	rand.seed(_init, 5);

	m_metallicity = system->m_metallicity * rand.Sfloat();
	m_population = sfloat(0);

	/* Bad type of planet for settlement */
	if (
		(averageTemp > CELSIUS+100) ||
		(averageTemp < 100) ||
	        ((type != SBody::TYPE_PLANET_DWARF) &&
		 (type != SBody::TYPE_PLANET_SMALL) &&
		 (type != SBody::TYPE_PLANET_WATER) &&
		 (type != SBody::TYPE_PLANET_CO2) &&
		 (type != SBody::TYPE_PLANET_METHANE) &&
		 (type != SBody::TYPE_PLANET_INDIGENOUS_LIFE) &&
		 (type != SBody::TYPE_PLANET_TERRAFORMED_POOR) &&
		 (type != SBody::TYPE_PLANET_TERRAFORMED_GOOD)
		)
	   )
	{
		return;
	}

	m_agricultural = sfloat(0);

	if ((type == SBody::TYPE_PLANET_INDIGENOUS_LIFE) ||
		(type == SBody::TYPE_PLANET_TERRAFORMED_GOOD)) {
		m_agricultural = CLAMP(sfloat(1,1) - sfloat(CELSIUS+25-averageTemp, 40), sfloat(0), sfloat(1,1));
		system->m_agricultural += 2*m_agricultural;
	} else if (type == SBody::TYPE_PLANET_TERRAFORMED_POOR) {
		m_agricultural = CLAMP(sfloat(1,1) - sfloat(CELSIUS+30-averageTemp, 50), sfloat(0), sfloat(1,1));
		system->m_agricultural += 1*m_agricultural;
	} else {
		// don't bother populating crap planets
		if (m_metallicity < sfloat(5,10)) return;
	}

	const int NUM_CONSUMABLES = 10;
	const Equip::Type consumables[NUM_CONSUMABLES] = { 
		Equip::AIR_PROCESSORS,
		Equip::GRAIN,
		Equip::FRUIT_AND_VEG,
		Equip::ANIMAL_MEAT,
		Equip::LIQUOR,
		Equip::CONSUMER_GOODS,
		Equip::MEDICINES,
		Equip::HAND_WEAPONS,
		Equip::NARCOTICS,
		Equip::LIQUID_OXYGEN
	};

	/* Commodities we produce (mining and agriculture) */
	for (int i=(int)Equip::FIRST_COMMODITY; i<(int)Equip::LAST_COMMODITY; i++) {
		Equip::Type t = (Equip::Type)i;
		const EquipType &type = EquipType::types[t];
		if (type.techLevel > system->m_techlevel) continue;

		sfloat affinity = sfloat(1,1);
		if (type.econType & ECON_AGRICULTURE) {
			affinity *= 2*m_agricultural;
		}
		if (type.econType & ECON_INDUSTRY) affinity *= system->m_industrial;
		// make industry after we see if agriculture and mining are viable
		if (type.econType & ECON_MINING) {
			affinity *= m_metallicity;
		}
		affinity *= rand.Sfloat();
		// producing consumables is wise
		for (int j=0; j<NUM_CONSUMABLES; j++) {
			if (i == consumables[j]) affinity *= 2; break;
		}
		assert(affinity >= 0);
		/* workforce... */
		m_population += affinity * system->m_humanProx;
		
		int howmuch = (affinity * 256).ToInt32();

		system->m_tradeLevel[t] += -2*howmuch;
		for (int i=0; i<EQUIP_INPUTS; i++) {
			if (!type.inputs[i]) continue;
			system->m_tradeLevel[type.inputs[i]] += howmuch;
		}
	}

	if (m_population > sfloat(1,10)) NameGenerator::PlanetName(rand);
	
	// Add a bunch of things people consume
	for (int i=0; i<NUM_CONSUMABLES; i++) {
		Equip::Type t = consumables[i];
		if ((t == Equip::AIR_PROCESSORS) ||
		    (t == Equip::LIQUID_OXYGEN) ||
		    (t == Equip::GRAIN) ||
		    (t == Equip::FRUIT_AND_VEG) ||
		    (t == Equip::ANIMAL_MEAT)) {
			if ((type == SBody::TYPE_PLANET_INDIGENOUS_LIFE) ||
			   (type == SBody::TYPE_PLANET_TERRAFORMED_GOOD)||
			   (type == SBody::TYPE_PLANET_TERRAFORMED_POOR))
				continue;
		}
		system->m_tradeLevel[t] += rand.Int32(32,128);
	}
	// well, outdoor worlds should have way more people
	m_population = sfloat(1,10)*m_population + m_population*m_agricultural;

//	printf("%s: pop %.3f billion\n", name.c_str(), m_population.ToFloat());

	outTotalPop += m_population;
}

void SBody::PopulateAddStations(StarSystem *system)
{
	for (unsigned int i=0; i<children.size(); i++) {
		children[i]->PopulateAddStations(system);
	}
	unsigned long _init[5] = { system->m_loc.systemNum, system->m_loc.sectorX,
			system->m_loc.sectorY, this->seed, UNIVERSE_SEED };
	MTRand rand;
	rand.seed(_init, 5);

	if (m_population < sfloat(1,1000)) return;

	sfloat pop = m_population + rand.Sfloat();

	sfloat orbMax = sfloat(1,4)*this->CalcHillRadius();
	sfloat orbMin = 4 * this->radius * AU_EARTH_RADIUS;
	if (children.size()) orbMax = MIN(orbMax, sfloat(1,2) * children[0]->orbMin);

	// starports - orbital
	pop -= rand.Sfloat();
	if ((orbMin < orbMax) && (pop >= 0)) {
	
		SBody *sp = system->NewBody();
		sp->type = SBody::TYPE_STARPORT_ORBITAL;
		sp->seed = rand.Int32();
		sp->tmp = 0;
		sp->parent = this;
		sp->rotationPeriod = sfloat(1,3600);
		sp->averageTemp = this->averageTemp;
		sp->mass = 0;
		sp->name = NameGenerator::Surname(rand) + " Spaceport";
		/* just always plonk starports in near orbit */
		sp->semiMajorAxis = orbMin;
		sp->eccentricity = sfloat(0);
		sp->axialTilt = sfloat(0);
		sp->orbit.eccentricity = 0;
		sp->orbit.semiMajorAxis = sp->semiMajorAxis.ToDouble()*AU;
		sp->orbit.period = calc_orbital_period(sp->orbit.semiMajorAxis, this->mass.ToDouble() * EARTH_MASS);
		sp->orbit.rotMatrix = matrix4x4d::Identity();
		children.insert(children.begin(), sp);
		system->m_spaceStations.push_back(sp);
		sp->orbMin = sp->semiMajorAxis;
		sp->orbMax = sp->semiMajorAxis;

		pop -= rand.Sfloat();
		if (pop > 0) {
			SBody *sp2 = system->NewBody();
			*sp2 = *sp;
			sp2->orbit.rotMatrix = matrix4x4d::RotateZMatrix(M_PI);
			sp2->name = NameGenerator::Surname(rand) + " Spaceport";
			children.insert(children.begin(), sp2);
			system->m_spaceStations.push_back(sp2);
		}
	}
	// starports - surface
	pop = m_population + rand.Sfloat();
	int max = 6;
	while (max-- > 0) {
		pop -= rand.Sfloat();
		if (pop < 0) break;

		SBody *sp = system->NewBody();
		sp->type = SBody::TYPE_STARPORT_SURFACE;
		sp->seed = rand.Int32();
		sp->tmp = 0;
		sp->parent = this;
		sp->averageTemp = this->averageTemp;
		sp->mass = 0;
		sp->name = NameGenerator::Surname(rand) + " Starport";
		memset(&sp->orbit, 0, sizeof(Orbit));
		position_settlement_on_planet(sp);
		children.insert(children.begin(), sp);
		system->m_spaceStations.push_back(sp);
	}
}

StarSystem::~StarSystem()
{
	if (rootBody) delete rootBody;
}

bool StarSystem::IsSystem(int sector_x, int sector_y, int system_idx)
{
	return (sector_x == m_loc.sectorX) && (sector_y == m_loc.sectorY) && (system_idx == m_loc.systemNum);
}

SBody::~SBody()
{
	for (std::vector<SBody*>::iterator i = children.begin(); i != children.end(); ++i) {
		delete (*i);
	}
}

void StarSystem::Serialize(Serializer::Writer &wr, StarSystem *s)
{
	if (s) {
		wr.Byte(1);
		wr.Int32(s->m_loc.sectorX);
		wr.Int32(s->m_loc.sectorY);
		wr.Int32(s->m_loc.systemNum);
	} else {
		wr.Byte(0);
	}
}

StarSystem *StarSystem::Unserialize(Serializer::Reader &rd)
{
	if (rd.Byte()) {
		int sec_x = rd.Int32();
		int sec_y = rd.Int32();
		int sys_idx = rd.Int32();
		return new StarSystem(sec_x, sec_y, sys_idx);
	} else {
		return 0;
	}
}

#define STARSYS_MAX_CACHED 8
static std::list<StarSystem*> s_cachedSystems;

StarSystem *StarSystem::GetCached(const SysLoc &loc)
{
	for (std::list<StarSystem*>::iterator i = s_cachedSystems.begin();
			i != s_cachedSystems.end(); ++i) {
		if ((*i)->m_loc == loc) {
			// move to front of cache to indicate it is hot
			StarSystem *s = *i;
			s_cachedSystems.erase(i);
			s_cachedSystems.push_front(s);
			return s;
		}
	}
	StarSystem *s = new StarSystem(loc.sectorX, loc.sectorY, loc.systemNum);
	s_cachedSystems.push_front(s);
	return s;
}

void StarSystem::ShrinkCache()
{
	int n=0;
	for (std::list<StarSystem*>::iterator i = s_cachedSystems.begin();
			i != s_cachedSystems.end(); ++i, n++) {
		if (n >= STARSYS_MAX_CACHED) {
			while (i != s_cachedSystems.end()) {
				delete *i;
				i = s_cachedSystems.erase(i);
			}
			break;
		}
	}
}
