#ifndef PRAD_ET_CHANNEL_H
#define PRAD_ET_CHANNEL_H

#include <unordered_map>
#include <memory>
#include <string>
#include <stdint.h>
#include "et.h"
#include "PRadException.h"
#include "PRadETStation.h"

#define ET_CHUNK_SIZE 500

class PRadETChannel
{
// configuration class
public:
    class Configuration
    {
    public:
        Configuration() {
            void *ptr;
            et_open_config_init(&ptr);
            conf = std::shared_ptr<void>(ptr, [] (void *p) { et_open_config_destroy(p); });
        }

        et_openconfig Get() { return conf.get(); }

        // wrapper functions
        void SetWait(int val) { et_open_config_setwait(conf.get(), val); }
        int GetWait() { int val; et_open_config_getwait(conf.get(), &val); return val; }

        void SetCast(int val) { et_open_config_setcast(conf.get(), val); }
        int GetCast() { int val; et_open_config_getcast(conf.get(), &val); return val; }

        void SetTTL(int val) { et_open_config_setTTL(conf.get(), val); }
        int GetTTL() { int val; et_open_config_getTTL(conf.get(), &val); return val; }

        void SetMode(int val) { et_open_config_setmode(conf.get(), val); }
        int GetMode() { int val; et_open_config_getmode(conf.get(), &val); return val; }

        void SetDebugDefault(int val) { et_open_config_setdebugdefault(conf.get(), val); }
        int GetDebugDefault() { int val; et_open_config_getdebugdefault(conf.get(), &val); return val; }

        void SetPort(int val) { et_open_config_setport(conf.get(), val); }
        int GetPort() { int val; et_open_config_getport(conf.get(), &val); return val; }

        void SetServerPort(int val) { et_open_config_setserverport(conf.get(), val); }
        int GetServerPort() { int val; et_open_config_getserverport(conf.get(), &val); return val; }

        void SetTimeOut(struct timespec val) { et_open_config_settimeout(conf.get(), val); }
        struct timespec GetTimeOut() { struct timespec val; et_open_config_gettimeout(conf.get(), &val); return val; }

        void SetHost(const std::string &val) { et_open_config_sethost(conf.get(), val.c_str()); }
        std::string GetHost() { char val[1024]; et_open_config_gethost(conf.get(), val); return std::string(val); }

        void AddBroadCast(const std::string &val) { et_open_config_addbroadcast(conf.get(), val.c_str()); }
        void RemoveBroadCast(const std::string &val) { et_open_config_removebroadcast(conf.get(), val.c_str()); }

        void AddMultiCast(const std::string &val) { et_open_config_addmulticast(conf.get(), val.c_str()); }
        void RemoveMultiCast(const std::string &val) { et_open_config_removemulticast(conf.get(), val.c_str()); }

        void SetPolicy(int val) { et_open_config_setpolicy(conf.get(), val); }
        int GetPolicy() { int val; et_open_config_getpolicy(conf.get(), &val); return val; }

        void SetInterface(const std::string &val) { et_open_config_setinterface(conf.get(), val.c_str()); }
        std::string GetInterface() { char val[1024]; et_open_config_getinterface(conf.get(), val); return std::string(val); }

        void SetTCP(int rBufSize, int sBufSize, int noDelay) { et_open_config_settcp(conf.get(), rBufSize, sBufSize, noDelay); }
        void GetTCP(int *rBufSize, int *sBufSize, int *noDelay) { et_open_config_gettcp(conf.get(), rBufSize, sBufSize, noDelay); }

    private:
        std::shared_ptr<void> conf;
    };

public:
    PRadETChannel(size_t size = 1048576);
    virtual ~PRadETChannel();

    void Open(const char *ipAddr, int tcpPort, const char *etFile);
    void NewStation(const std::string &name, int mode=2);
    void SwitchStation(const std::string &name);
    void RemoveStation(const std::string &name);
    void AttachStation();
    void DetachStation();
    void ForceClose();
    bool Read();
    bool Write(void *buf, int nbytes);
    void *GetBuffer() {return (void*) buffer;}
    size_t GetBufferLength() {return bufferSize;}
    Configuration &GetConfig() {return config;}
    et_sys_id &GetID() {return et_id;}
    PRadETStation *GetCurrentStation() {return curr_stat;}
    PRadETStation *GetStation(const std::string &name);

private:
    Configuration config;
    PRadETStation *curr_stat;
    std::unordered_map<std::string, PRadETStation*> stations;
    et_sys_id et_id;
    et_event *etEvent;
    uint32_t *buffer;
    size_t bufferSize;
    void copyEvent();

public:
    static std::string get_error_str(int error)
    {
        switch(error) {
        case ET_ERROR: return "General error.";
        case ET_ERROR_TOOMANY: return "Too many somethings (stations, attachments, temp events, ET system responses) exist.";
        case ET_ERROR_EXISTS: return "ET system file or station already exists.";
        case ET_ERROR_WAKEUP: return "Sleeping routine woken up by et_wakeup_attachment() or et_wakeup_all().";
        case ET_ERROR_TIMEOUT: return "Timed out.";
        case ET_ERROR_EMPTY: return "No events available in async mode.";
        case ET_ERROR_BUSY: return "Resource is busy.";
        case ET_ERROR_DEAD: return "ET system is dead.";
        case ET_ERROR_READ: return "Network read error.";
        case ET_ERROR_WRITE: return "Network write error,";
        case ET_ERROR_REMOTE: return "Cannot allocate memory in remote client.";
        case ET_ERROR_TOOBIG: return "Client is 32 bits & server is 64 (or vice versa) and event is too big for one.";
        case ET_ERROR_NOMEM: return "Cannot allocate memory.";
        case ET_ERROR_BADARG: return "Bad argument given to function.";
        case ET_ERROR_SOCKET: return "Socket option could not be set.";
        case ET_ERROR_NETWORK: return "Host name or address could not be resolved, or cannot connect.";
        case ET_ERROR_CLOSED: return "ET system has been closed by client.";
        case ET_ERROR_JAVASYS: return "C code trying to open Java-based ET system file locally.";
        default: break;
        }
        return "Unknown error";
    }
};

#endif
